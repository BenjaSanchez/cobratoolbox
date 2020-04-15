function [annotationString,notes] = makeSBMLAnnotationString(model,id,fieldentries,position,cobrapy)
% makeSBMLAnnotationString gives the annotationString for an SBML based on the fields in the model
%
% USAGE:
%
%       [annotationString, notes] = makeSBMLAnnotationString(model,id,fieldentries,position)
%
% INPUT:
%    model:            the model to extract the data
%    id:               the ID of the entity
%    fieldentries:     either a char indicating the field
%                      (prot,met,rxn,comp,gene), or a cell array with X{:,1}
%                      being field IDs and X{:,2} being bioql qualiiers to
%                      annotate for the field.
%    position:         the position in the model to extract the data.
%    cobrapy:          Boolean to guarantee compatibility with cobrapy (default=false)
%
% OUTPUT:
%
%   annotationString: The annotation String to be put into the SBML.
%   notes:            A 2*x cell array of fields which did not contain
%                     valid identifiers (according to the pattern check.
%
% .. Authors:
%       - Thomas Pfau May 2017 



if nargin < 5
    cobrapy = false;
end

allQualifiers = getBioQualifiers();

if cobrapy
    allQualifiers(strcmp(allQualifiers,'hasProperty')) = [];  %cobrapy does not store SBO terms as annotations (already included)
end
if ischar(fieldentries)
    fieldentries = {fieldentries, allQualifiers};
end
annotationString = '';
modelFields = fieldnames(model);
%Only look at the relevant fields.
tmp_note = '';
notes = cell(0,2);
bagindentlevel = '      ';
for pos = 1:size(fieldentries,1)        
    field = fieldentries{pos,1};
    if isempty(fieldentries{pos,2})
        allowedQualifiers = allQualifiers;
    else
        allowedQualifiers = fieldentries{pos,2};
    end
    fieldMappings = getDatabaseMappings(field);
    [~,upos,~] = unique(fieldMappings(:,3));
    fieldMappings = fieldMappings(upos,:);
    
    relfields = modelFields(cellfun(@(x) strncmp(x,field,length(field)),modelFields));
    
    for i = 1:numel(allowedQualifiers)
        annotationsFields = relfields(cellfun(@(x) strncmp(x,[field allowedQualifiers{i}],length([field allowedQualifiers{i}])),relfields));
        knownFields = fieldMappings(cellfun(@(x) strcmp(x,allowedQualifiers{i}),fieldMappings(:,2)),:);
        dbnote = '';
        dblist = {};
        for fieldid = 1:numel(annotationsFields)
            if isempty(model.(annotationsFields{fieldid}){position})
                continue
            end
            ids = strsplit(model.(annotationsFields{fieldid}){position},';');
            
            
            dbname = convertSBMLID(regexprep(annotationsFields{fieldid},[field allowedQualifiers{i} '(.*)' 'ID$'],'$1'),false);
            dbrdfstring = [bagindentlevel '    <rdf:li rdf:resource="https://identifiers.org/' dbname '/'];
            if cobrapy
                ids = strrep(ids, ' ', '');  % bug in some ids
            end
            dbstring = strjoin(strcat(dbrdfstring,ids,sprintf('%s\n','"/>')),sprintf('\n'));
            if cobrapy
                dblist = [dblist; dbstring];
            else
                dbnote = [dbnote, dbstring];
            end
        end
        knownExistentFields = knownFields(ismember(knownFields(:,3),modelFields),:);
        
        for fieldid = 1:size(knownExistentFields,1)
            if isempty(model.(knownExistentFields{fieldid,3}){position})
                continue
            end
            ids = strtrim(strsplit(model.(knownExistentFields{fieldid,3}){position},';'));
            correctids = ~cellfun(@isempty, regexp(ids,knownExistentFields{fieldid,5}));            
            %If we have correct ids, we will annotate those.
            if any(correctids)
                dbname = knownExistentFields{fieldid,1};
                dbrdfstring = [bagindentlevel '    <rdf:li rdf:resource="https://identifiers.org/' dbname '/'];
                dbstring = strjoin(strcat(dbrdfstring,ids(correctids),sprintf('%s\n','"/>')),sprintf('\n'));
                if cobrapy
                    dblist = [dblist; dbstring];
                else
                    dbnote = [dbnote, dbstring];
                end
            end
            %if we have incorrect ids, we will add this data to the notes
            %of the reaction.
            if any(~correctids)
                notes(end+1,:) = {knownExistentFields{fieldid,3}, model.(knownExistentFields{fieldid,3}){position}};
            end
        end
        
        if cobrapy
            annotationsFields = regexprep(annotationsFields,'^metis','');
            annotationsFields = regexprep(annotationsFields,'^rxnis','');
            annotationsFields = regexprep(annotationsFields,'ID$','');
            annotationsFields = convertSBMLID(annotationsFields, false);
            allFields = sort([annotationsFields; knownExistentFields(:,1)]);
            for fieldid = 1:numel(allFields)
                for dbid = 1:numel(dblist)
                    if contains(dblist{dbid}, allFields{fieldid})
                        dbnote = [dbnote, dblist{dbid}];
                    end
                end
            end
        end
        if ~isempty(dbnote)
            %Make specification for this bag
            specstring = ['<bqbiol:' allowedQualifiers{i} ' xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">' ];
            specend = ['</bqbiol:' allowedQualifiers{i} '>'];
            
            specification_string = sprintf('%s%s\n%s  %s\n%s%s\n',bagindentlevel,specstring,bagindentlevel,'<rdf:Bag/>',bagindentlevel,specend);                        
             
            tmp_note=[specification_string, tmp_note, bagindentlevel, '<bqbiol:', allowedQualifiers{i}, sprintf('%s\n%s  %s\n','>', bagindentlevel,'<rdf:Bag>')];
            tmp_note = [tmp_note, dbnote ];
            tmp_note = [ tmp_note, sprintf('\n  %s%s\n%s%s\n',bagindentlevel,'</rdf:Bag>',bagindentlevel, ['</bqbiol:' allowedQualifiers{i} '>'])]; % ending syntax
        end
        
    end
end
if ~isempty(tmp_note)
    if cobrapy
        annotopentag = '<annotation>';
    else
        annotopentag = '<annotation xmlns:sbml="http://www.sbml.org/sbml/level3/version1/core">';
    end
    rdfOpenTag = '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:vCard4="http://www.w3.org/2006/vcard/ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">';
    annotationString = sprintf('%s\n  %s\n    ',annotopentag,rdfOpenTag);   
    annotationString = [ annotationString '<rdf:Description rdf:about="#',id,'">'];
    annotationString = [annotationString sprintf('\n')];
    annotationString = [annotationString, tmp_note, sprintf('    %s\n  %s\n%s','</rdf:Description>', '</rdf:RDF>', '</annotation>')];
end
end
