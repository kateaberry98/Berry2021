function  Neighbours = Moore(J,X,Y)
% Records indices of nearest neighbours for products/shelves

Neighbours = [J+Y,J+Y+1,J+1,J-Y+1,J-Y,J-Y-1,J-1,J+Y-1];
Neighbours = Neighbours.*(Neighbours > 0 & Neighbours < (X*Y+1));
Neighbours = Neighbours(Neighbours > 0);

end