function [Xtrans,Ytrans,Angle] = registration(imgTriangle)
  middleLine = size(imgTriangle,2)/2;
  imgRidges  = imgTriangle(:,1:middleLine);
  imgContact = imgTriangle(:,middleLine+1:end);

  imgContactCorr = imgContact;
  imgContactCorr(imgContactCorr <0) = 0;

  imgRidgesCorr = imgRidges;
  imgRidgesCorr(imgRidgesCorr < 0) = 0;

  t=[5 5];
  H2 = fspecial('average',t);
  ImageMoy = imfilter(imgContactCorr, H2); 
  ImageMoy = medfilt2(imgContactCorr, t);

  imgContactCorr = imadjust(imgContactCorr);
  imgRidgesCorr = imadjust(imgRidgesCorr);

  [centers1,radii1] = imfindcircles(ImageMoy,[20 40],'ObjectPolarity','dark', ...
      'Sensitivity',0.925);
  %     figure(1)
  %     imshow(ImageMoy)
  %     h = viscircles(centers1,radii1);

  [centers2,radii2] = imfindcircles(imgRidgesCorr,[20 40],'ObjectPolarity','dark', ...
      'Sensitivity',0.95);
  %   figure(1); hold on
  %   imshow([imgRidgesCorr,ImageMoy])
  %   h=viscircles(centers2,radii2);
  %   h=viscircles(centers1+[middleLine,0],radii1);

  centers1XSorted = sort (centers1(:,1));
  centers1YSorted = sort (centers1 (:,2));

  centers2XSorted = sort (centers2(:,1));
  centers2YSorted = sort (centers2 (:,2));
  Xtrans =  -(centers1XSorted(3,1)-centers2XSorted(3,1)+7);
  Ytrans =  -(centers1YSorted(2,1)-centers2YSorted(2,1));

  tan1 = (centers1YSorted(2,1)-centers1YSorted(1,1)) / (centers1XSorted(2,1)-centers1XSorted(1,1)) ;
  tan2 = (centers2YSorted(2,1)-centers2YSorted(1,1)) / (centers2XSorted(2,1)-centers2XSorted(1,1)) ;
  atan1 = (atan (tan1));
  atan2 = (atan (tan2)) ;

  Angle = (atan1 - atan2) ;
  Angle = ( Angle * 180 ) / pi ;
end
