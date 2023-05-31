%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Load mire and blank image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%compute scale with one coin
imgPiece = im2double(imread([folder img_folder '/piece.jpg']));
scale = scaleComputation(imgPiece); % in mm/px

%%estimate the transformation and realign figures
imgTriangle = im2double(imread([folder img_folder '/triangle1.jpg']));
[Xtrans,Ytrans,Angle] = registration(imgTriangle);

middleLine = size(imgTriangle,2)/2;
imgRidges = imgTriangle(:,1:middleLine);
imgContact = imgTriangle(:,middleLine+1:end);
imgContactRot = imrotate (imgContact, Angle , 'crop');
imgContactReg = imtranslate(imgContactRot,[Xtrans ,Ytrans]);

%%load background and correct contrast and light deviation
imgBackground = im2double(imread([folder img_folder '/vide.jpg']));
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

