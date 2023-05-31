function scale = scaleComputation(imgPiece)
  middleLine = size(imgPiece,2)/2;
  imgContact = imgPiece(:,middleLine+1:end);
  [centers1,radii1] = imfindcircles(imgContact,[50 100],'ObjectPolarity','dark', ...
    'Sensitivity',0.975);
  scale = 2.75/radii1; % radius of the real piece = 2.75
end
