function [pathPairs ] = convertFromCuttingToMesh( Segs )
%CONVERTFROMCUTTINGTOMESH Summary of this function goes here
%   Detailed explanation goes here
M=TriangleMesh('VF',Segs.V,Segs.F);
E=M.E;
% B=M.findOrientedBoundaries;
% Bm=zeros(numel(B),max(cellfun(@numel,B)));
% for i=1:numel(B)
%     Bm(i,1:numel(B{i}))=B{i};
% end
O2N=Segs.Old2New;
Path=Segs.Paths;
pathPairs={};%cell(size(Path));
S=Segs.Segments;
for i=1:size(Path,1)
    for j=1:size(Path,2)
        if isempty(Path{i,j})
            continue;
        end
        Pi=O2N(Path{i,j});
        edges=sort([Pi{1};repmat(Pi{2}(1),1,numel(Pi{1}))]',2);
        Pi1=Pi{1}(find(ismember(edges,E','rows')));
        
        edges=sort([Pi{1};repmat(Pi{2}(2),1,numel(Pi{1}))]',2);
        Pi2=Pi{1}(find(ismember(edges,E','rows')));
        Pi{1}=[Pi1,Pi2];
        
        edges=sort([Pi{end};repmat(Pi{end-1}(1),1,numel(Pi{end}))]',2);
        Pi1=Pi{end}(find(ismember(edges,E','rows')));
        
        edges=sort([Pi{end};repmat(Pi{end-1}(2),1,numel(Pi{end}))]',2);
        Pi2=Pi{end}(find(ismember(edges,E','rows')));
        Pi{end}=[Pi1,Pi2];
        
        %     C1=S(Path1icut{17}(1));
        %     C2=S(Path1icut{17}(2));
        %     ind1=find(S(Path1icut{1})==C1);
        %
        %     ind2=find(S(Path1icut{1})==C2);
        %     Path1icut{1}=Path1icut{1}([ind1,ind2]);
        %
        %     ind1=find(S(Path1icut{end})==C1);
        %     ind2=find(S(Path1icut{end})==C2);
        %     Path1icut{end}=Path1icut{end}([ind1,ind2]);
%         pathPairs{i,j}=cell2mat(Pi');
        pathPairs{end+1}=cell2mat(Pi');
    end
end
