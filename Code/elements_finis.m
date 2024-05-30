function elements_finis(r1,r2,V1,V2,a,f,Hmax)
if nargin == 0
    r1=0.2;r2=1;V1=12;V2=10;a=10;f=-10;Hmax=0.01;
end

close all

% Create geometry model
model = createpde;

C1 = [1 0 0 r1]';C2 = [1 0 0 r2]';
geom = [C1 C2];

ns = (char('C1','C2'))';

% Set formula
sf = 'C2 - C1';

% Create geometry
g = decsg(geom,sf,ns);
% g(2,:).^2+g(4,:).^2
% g(3,:).^2+g(5,:).^2

geometryFromEdges(model,g);

figure,pdegplot(model,'EdgeLabels','on'),drawnow

specifyCoefficients(model,'m',0,'d',0,'c',1,'a',a,'f',f);

applyBoundaryCondition(model,'dirichlet','edge',[1 2 3 4],'u',V1);
applyBoundaryCondition(model,'dirichlet','edge',[5 6 7 8],'u',V2);

mesh = generateMesh(model,'Hmax',Hmax,'GeometricOrder','linear');

results = solvepde(model);
u_n = results.NodalSolution;
disp(length(u_n))

figure,pdeplot(model)
drawnow

figure,pdeplot(model,'XYData',u_n,'ZData',u_n)
colormap(jet)
drawnow

labels=zeros(size(u_n));
labels(u_n==V1)=1;
labels(u_n==V2)=2;

params=[r1 r2 V1 V2 a f];

% save elecdata
writetable(table(mesh.Nodes'),'nodes.csv')
writetable(table(mesh.Elements'),'elements.csv')
writetable(table(u_n),'solution.csv')
writetable(table(labels),'labels.csv')
writetable(table(params),'params.csv')