function map=FlatteningsToMap(fs,ft)
V_s=fs.flat_V;
V_t=ft.flat_V;
T_s=fs.flat_T;
T_t=ft.flat_T;
pathPairs_s=fs.M_cut.pathPairs;
pathPairs_t=ft.M_cut.pathPairs;
gS={};
T={};
pathEnds=[];
for i=1:length(pathPairs_s)
    p1=V_s(pathPairs_s{i}(1,1),:);
    p2=V_s(pathPairs_s{i}(end,1),:);
    q1=V_s(pathPairs_s{i}(1,2),:);
    q2=V_s(pathPairs_s{i}(end,2),:);
    [A,t]=computeSimilarity(p1,p2,q1,q2);
    gS{i}=A;
    T{i}=t';
    pathEnds=[pathEnds;reshape(pathPairs_s{i}([1 end],:),4,1) reshape(pathPairs_t{i}([1 end],:),4,1)];
end
pathEnds=unique(pathEnds,'rows');
source=[];
source.TR=triangulation(T_s,V_s);
source.Y=V_s;
source.pathPairs=pathPairs_s;
target=[];
target.TR=triangulation(T_t,V_t);
target.Y=V_t;
target.pathPairs=pathPairs_t;
lifter1to2=RobustLifter(source,target,gS,T,pathEnds);
lifter1to2.lift();
lifter2to1=RobustLifter(target,source,gS,T,pathEnds(:,[2 1]));
lifter2to1.lift();
map=UncutSurfMap(fs.M_orig,ft.M_orig,fs.M_cut,ft.M_cut,lifter1to2,lifter2to1);
end