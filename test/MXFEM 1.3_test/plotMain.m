% Written By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function plotMain(omega,DISP)
% This function defines which plots are to be made based on user-defined
% PLOT matrix.

global DOMAIN PLOT

for iPlot = 1:5
    if PLOT(iPlot,1) == 1
        if iPlot == 1                                                       % Plot level set
            plotLevelSet;
        elseif iPlot == 2                                                   % Plot mesh
            plotMesh;
        elseif iPlot == 3                                                   % Plot deformation
            nXElem = DOMAIN(1);                                             % Number of elements in x-direction
            nYElem = DOMAIN(2);                                             % Number of elements in y-direction
            plotDeformation(DISP(1:(2*(nXElem+1)*(nYElem+1))));
        elseif iPlot == 4                                                   % Plot stress
            [Sxx,Sxy,Syy,Svm] = elemStress(omega,DISP);
            plotStress(Sxx,Sxy,Syy,Svm);
        elseif iPlot == 5                                                   % Plot stress contour
            if PLOT(4,1) ~= 1
                [Sxx,Sxy,Syy,Svm] = elemStress(omega,DISP);
            end
            plotContour(Sxx,Sxy,Syy,Svm);
        end
    end
end