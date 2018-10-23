rhoW = 19600; rhoE = 1060;
fracW = 0.87; fracE = 1-fracW;
rhoT = (rhoW*rhoE)/(fracW*rhoE+fracE*rhoW);
disp(['Total density is' rhoT 'kg/m3']) 