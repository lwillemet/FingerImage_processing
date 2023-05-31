%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Load mire and blank image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = pwd;
%%compute scale with one coin
imgPiece = im2double(imread([folder '/piece.jpg']));
scale = scaleComputation(imgPiece); % in mm/px

%%estimate the transformation and realign figures
imgTriangle = im2double(imread([folder '/triangle1.jpg']));
[Xtrans,Ytrans,Angle] = registration(imgTriangle);

middleLine = size(imgTriangle1,2)/2;
imgRidges = imgTriangle(:,1:middleLine);
imgContact = imgTriangle(:,middleLine+1:end);
imgContactRot = imrotate (imgContact, Angle , 'crop');
imgContactReg = imtranslate(imgContactRot,[Xtrans ,Ytrans]);

%%load background and correct contrast and light deviation
imgBackground = im2double(imread([folder '/empty.jpg']));
imgRidgesBack  = imgBackground(:,1:middleLine);
imgContactBack  = imgBackground(:,middleLine+1:end);
imgContactRot = imrotate (imgContactBack, Angle , 'crop');
imgContactBackReg = imtranslate(imgContactRot,[Xtrans ,Ytrans]);
histogram(imgContactBackReg(50:450,150:600),0:0.003:0.5)


    %%transmissive light, black is background
imgRidgesCorr = imgRidges-imgRidgesBack;
imgRidgesCorr(imgRidgesCorr <0) = 0;
imgRidgesCorr = imadjust(imgRidgesCorr,[0,0.5]);

    %%reflective light, white is background
    %%scale according to background illumination
se = strel('disk',9);
% imgContactBackOpen = (imopen(imgContactBackReg,se)+0.15)*1.15;
imgContactBackOpen = imgContactBackReg;
imgContactCorr = imgContactReg./imgContactBackOpen;
%imgContactCorr = imgContactCorr/prctile(imgContactCorr(:),99);  %%scale back to 1
imgContactCorr(imgContactCorr >1) = 1;                          %%clip extreme values

%%find the circle in which the circle enters in contact
h = imshow(imgContactReg);
e = imellipse(gca,[128 54 482 420]); % must be adjusted for each participant
imgMaskContact = createMask(e,h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoImg = imfinfo(filenameImg);

fsImg = 1000; %%one thousand frames per second
nb_img_before_trigger = 300;
t1 = nb_img_before_trigger/fsImg; %time instant of the first image
tImg = -t1:1/fsImg:(length(infoImg)-nb_img_before_trigger-1)/fsImg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Load and process force data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(filenameData)

[idxTrig] = find(trigger>mean(trigger),1,'first');
if idxTrig <=10
    [~,idxTrig] = min(trigger);
end
tData = time-time(idxTrig);
fsData = 31e3;

forceInterp = abs(interp1(tData,rawForces(3,:),tImg,'linear'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Load and process images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff = 0.01;
numFrameTrack = find(forceInterp>ff,1,'first'); %%find the reference points when skin is present

lastFrame = length(tImg);

%%load and filter all pictures
for kk = numFrameTrack:lastFrame
    
    %%read current img
    img = im2double(imread(filenameImg, kk, 'Info', infoImg));
    
    %%cut images and align
    middleLine = size(img,2)/2;
    imgRidges  = img(:,1:middleLine);
    imgContact = img(:,middleLine+1:end);
    imgContactRot = imrotate (imgContact, Angle , 'crop');
    imgContactReg = imtranslate(imgContactRot,[Xtrans ,Ytrans]);
    
    %%background removal
    imgRidgesCorr = imgRidges - imgRidgesBack;
    
    %%Image filtering
    imgRidgesAll(:,:,kk) = imadjust(imgRidgesCorr,[0 0.5]);
    imgSharp(:,:,kk) = imadjust(imsharpen(imgRidgesAll(:,:,kk),'Radius',5,'Amount',2),[.15,1]);   %%for illustration purposes

    imgContactRaw(:,:,kk) = imgContactReg./imgContactBackOpen;

end
[~,iFmax] = max(forceInterp);

for kk=numFrameTrack:lastFrame
    %%subtraction with first image
    imgContact = imgContactRaw(:,:,numFrameTrack).*imgMaskContact - imgContactRaw(:,:,kk).*imgMaskContact;
    imgContact(imgContact<0) = 0;
    imgContactAll(:,:,kk) = medfilt2(im2double(imgContact),[5 5]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute the real and the gross contact area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seD = strel('diamond',1);
se = strel('disk',8);
[~,idx1] = find(forceInterp>3,1,'first');
if (isempty(idx1))
    idx1 = iFmax;
end
thrC = 0.95*graythresh(imgContactAll(:,:,idx1));
thrMax = graythresh(imgContactAll(:,:,iFmax));
for kk = numFrameTrack:lastFrame
    
    %%real area of contact: made by the micro-junctions
    IA1 = imbinarize(imgContactAll(:,:,kk),thrC);
    imgContactBW(:,:,kk) = imerode(IA1,seD); %erode image with diamond
    
    metricBrightness(kk) = nansum(nansum(imgContactAll(:,:,kk)));
    metricArea(kk) = nansum(nansum(imgContactBW(:,:,kk)))*scale.^2;
    
    %%gross area of contact (ellipse shape)
    IA2 = imbinarize(imgContactAll(:,:,kk),thrMax);
    BWdil = imdilate(IA2,se);
    BWfill(:,:,kk) = imfill(BWdil,'holes');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Points tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extract ellipse of contact to track points of interest
imshow(BWfill(:,:,iFmax))
contact_ellipse = regionprops('table',BWfill(:,:,iFmax),'Centroid','MajorAxisLength','MinorAxisLength');

[~,ell_ii] = max([contact_ellipse.MajorAxisLength]);
xelipse=contact_ellipse.Centroid(ell_ii,1);
yelipse=contact_ellipse.Centroid(ell_ii,2);
grandaxeelipe=contact_ellipse.MajorAxisLength(ell_ii);
petitaxeelipe=contact_ellipse.MinorAxisLength(ell_ii);
h = imshow(BWfill(:,:,iFmax));
axialROI = images.roi.Ellipse(gca,'Center',[xelipse yelipse],...
    'Semiaxes',[grandaxeelipe/2 petitaxeelipe/2]);
axialROImask=createMask(axialROI,h);
close;

nb_points = 1200; % number of points to track
pointsAll = NaN(nb_points,2,lastFrame);

kk=lastFrame;
imgTarget = imgRidgesAll(:,:,kk).*axialROImask;
point = detectMinEigenFeatures(imgTarget,'MinQuality',1e-6);

tracker = vision.PointTracker('MaxBidirectionalError',2,'BlockSize',[101 101]);
initialize(tracker,point.selectStrongest(nb_points).Location,imgTarget);

[points,validity] = tracker.step(imgTarget);
if (validity(validity==0))
    disp('some shit is going on')
end
clear pointsAll brightnessPts
pointsAll(:,:,lastFrame) = points;

%%initialize displacement vectors
displacementpx(:,:,1) = [0,0];
vector(:,:,kk) = [0*ones(50,1),0*ones(50,1)];

h = fspecial('average',[50 50]);
        
for kk = lastFrame-1:-1:numFrameTrack
    imgTarget = imgRidgesAll(:,:,kk).*axialROImask;

    [points,validity] = tracker.step(imgTarget);
    if (validity(validity==0))
        disp('some shit is going on')
    end
    pointsAll(:,:,kk) = points;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Displacements computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Consider points once they are in contact
contactIdx = lastFrame.*ones(1,length(pointsAll));
brightness = zeros(nb_points,lastFrame);
dispPoints = zeros(nb_points,2,lastFrame);

for kk=numFrameTrack:lastFrame
    idx = find(~isnan(pointsAll(:,1,kk)));
    imgblured(:,:,kk) = imadjust(imgaussfilt(imgContactAll(:,:,kk),5),[0 0.1]);

    for ii=1:length(idx)
        if(round(pointsAll(idx(ii),2,kk))<size(BWfill,1) & round(pointsAll(idx(ii),1,kk))<size(BWfill,2))
        brightness(idx(ii),kk) = BWfill(round(pointsAll(idx(ii),2,kk)),round(pointsAll(idx(ii),1,kk)),kk);
        else
        brightness(idx(ii),kk) = NaN;
        end
    end
end
for ii=1:nb_points
if(any(brightness(ii,:)>0))
    contactIdx(ii) = find(brightness(ii,:)>mean(brightness(ii,:)),1,'first');
    dispPoints(ii,:,:) = pointsAll(ii,:,:) - pointsAll(ii,:,contactIdx(ii));
    dispPoints(ii,1,:) = smooth(dispPoints(ii,1,:),20);
    dispPoints(ii,2,:) = smooth(dispPoints(ii,2,:),20);
    dispPoints(ii,:,numFrameTrack:contactIdx(ii)-1) = zeros(size(dispPoints(ii,:,numFrameTrack:contactIdx(ii)-1)));
else
    dispPoints(ii,:,:) = zeros(1,2,lastFrame);
end
end

for kk=numFrameTrack:lastFrame
dispTotal(kk,1) = sqrt(nanmedian(dispPoints(:,1,kk)).^2 + nanmedian(dispPoints(:,1,kk)).^2).*scale;
end
dispTotal(1:numFrameTrack-1) = NaN(1,numFrameTrack-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Divergence and curl of the deformation field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xsq,Ysq] = meshgrid(1:size(imgRidgesAll(:,:,lastFrame),2),1:size(imgRidgesAll(:,:,lastFrame),1));
Xq = Xsq.*axialROImask;
Yq = Ysq.*axialROImask;

for kk=numFrameTrack+1:size(pointsAll,3)
    
    div(:,kk) = gradient(dispPoints(~isnan(pointsAll(:,1,kk)),1,kk),pointsAll(~isnan(pointsAll(:,1,kk)),1,kk))+...
        gradient(dispPoints(~isnan(pointsAll(:,1,kk)),2,kk),pointsAll(~isnan(pointsAll(:,1,kk)),2,kk));
    
    % Interpolation on this grid of points
    X = double(pointsAll(:,1,kk));
    Y = double(pointsAll(:,2,kk));
    div_clean = rmoutliers([X,Y,double(div(:,kk))]);
    div_interp(:,:,kk) = griddata(div_clean(:,1),div_clean(:,2),div_clean(:,3),Xq,Yq,'cubic');
    Id = ones(size(BWfill(:,:,kk)));
    Id(BWfill(:,:,kk)==0) = NaN;
    div_interp_mask(:,:,kk) = div_interp(:,:,kk).*Id;
    div_median(:,kk) = nanmedian(div_interp_mask(:,:,kk),'all');
    
    % Curl
    curl(:,kk) = gradient(dispPoints(~isnan(pointsAll(:,2,kk)),1,kk),pointsAll(~isnan(pointsAll(:,1,kk)),1,kk))-...
        gradient(dispPoints(~isnan(pointsAll(:,1,kk)),1,kk),pointsAll(~isnan(pointsAll(:,2,kk)),2,kk));
    curl_clean = rmoutliers([X,Y,double(curl(:,kk))]);
    curl_interp(:,:,kk) = griddata(curl_clean(:,1),curl_clean(:,2),curl_clean(:,3),Xq,Yq,'cubic');
    curl_interp_mask(:,:,kk) = curl_interp(:,:,kk).*Id;
    curl_median(:,kk) = nanmedian(curl_interp_mask(:,:,kk),'all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Strain computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Building the Delaunay triangulation
triD = delaunayTriangulation(double(pointsAll(:,:,lastFrame)));

for kk = numFrameTrack:lastFrame
    
    triD.Points = double(pointsAll(:,:,kk));
    trisize = triD.size(1);
    
    for pp = 1:trisize
        
        uv = double(pointsAll(triD.ConnectivityList(pp,:),:,kk)-...
            pointsAll(triD.ConnectivityList(pp,:),:,numFrameTrack));
        xy = double(pointsAll(triD.ConnectivityList(pp,:),:,kk));
        
        uMat = [uv(2,1)-uv(1,1);uv(3,1)-uv(2,1)];
        vMat = [uv(2,2)-uv(1,2);uv(3,2)-uv(2,2)];
        xyMat = [[xy(2,1)-xy(1,1),xy(2,2)-xy(1,2)];...
            [xy(3,1)-xy(2,1),xy(3,2)-xy(2,2)]];
        
        [du] = inv(xyMat)*uMat;
        [dv] = inv(xyMat)*vMat;
        
        epsxx(pp,kk) = du(1)+1/2*(du(1).^2+dv(1).^2);
        epsyy(pp,kk) = dv(1)+1/2*(du(2).^2+dv(2).^2);
        epsxy(pp,kk) = 1/2*(du(2)+dv(1)) + 1/2*(du(1)*du(2)+dv(1)*dv(2));
        
        centerPts(pp,:,kk) = triD.incenter([pp]);

    end
    strain(1:trisize,1:2,kk) = centerPts(1:trisize,1:2,kk);
    strain(1:trisize,4,kk) = epsxx(1:trisize,kk);
    strain(1:trisize,5,kk) = epsyy(1:trisize,kk);
    strain(1:trisize,6,kk) = epsxy(1:trisize,kk);

end
