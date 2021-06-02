% **********************Garrett Simpson**********************************
%         
%                       Version 10.05.2020-Final
%
% Input will be DICOM image series accompanied by the RT structure set.
% Program will output original data range as well as the dynamic limited
% range as defined by Collewet et al.                                                                       
%
% Program calculates the texture features using Ng = 64 with 
% Histogram equalization 
% Each image will have texture features calculated and output.
% 
% 
%
% Remove Isocenter structures or else error will be thrown 
% 
% 
%**************************************************************************

%***Begin Code***
tic
basedir = 'C:\Users\gns24\Box Sync\Remote Desktop\Desktop\PancPtData\Panc36';
files = dir(basedir);
fractions = files(3:size(files,1),1);
allp = {'','','','',''};
for fxNum = 1:size(fractions,1)
    test = dir(strcat(basedir,'\',fractions(fxNum).name));
  for q1 = 1:size(test,1)
        if contains(test(q1).name,'Images')
            imagepath = strcat(basedir,'\',fractions(fxNum).name,'\',test(q1).name);
        end
        if contains(test(q1).name,'RTst')
            structpath = strcat(basedir,'\',fractions(fxNum).name,'\',test(q1).name);
        end
  end
 
  tempImageFilePaths = dir(imagepath);
  dicomImageFileNames = "";
  for q2 =1:size(dir(imagepath),1)
     if (contains(tempImageFilePaths(q2).name,'2.16' )|| contains(tempImageFilePaths(q2).name,'MRI2.16' ))
         dicomImageFileNames = [dicomImageFileNames; tempImageFilePaths(q2).name]; 
     end
  end
  % Loading DICOM information
  img = LoadDICOMImages(imagepath,dicomImageFileNames);
  % Getting RTst file
  rtstfileName = dir(structpath);
  x = LoadDICOMStructures(structpath,rtstfileName(3).name,img);
  for z = 1:length(x)
    if(1 == (strcmp(string(x{z}.name),"GTVp")) ||1 == (strcmp(string(x{z}.name),"GTV")) ||1 == (strcmp(string(x{z}.name),"GTVm"))||1 == (strcmp(string(x{z}.name),"GTV_VR")))
        tempstruct{1,1} = struct(x{z});
        tempstruct{1,1}.name
        break
    end
  end
  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!! Check image resolutions !!!!!!!!!!!!!!!!!
volume = img.data;
scanType = 'MRscan';
pixelW = img.width(1)*10;
sliceS = img.width(3)*10;
R = 1;
scale = 'pixelW';
textType = 'Matrix';
quantAlgo = {'Equal'};
Ng = cell2mat({64});
xlfilename = 'C:\Users\gns24\Box Sync\Remote Desktop\Desktop\PancUpdatePatients.xlsx';
positions = ['C2:J2'];
positionstwo = ['K2:M2'];
positionsthree = ['N2:Z2'];
positionsfour = ['AA2:AM2'];
positionsfive = ['AN2:AR2'];
mask = tempstruct{1}.mask;
[ROIonly,levels] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,textType,'Equal',64);  
[GLCM] = getGLCM(ROIonly,levels);
[textures] = getGLCMtextures(getGLCM(ROIonly,levels));
[globaltexts] = getGlobalTextures(ROIonly,max(levels(:)));
[GLRLMtextures] = getGLRLMtexturesUpdate(getGLRLM(ROIonly,levels),round((tempstruct{1,1}.volume)/(img.width(1)*img.width(1)*img.width(3))));
[GLSZMtextures] = getGLSZMtextures(getGLSZM(ROIonly,levels));
[NGTDM,countValid] = getNGTDM(ROIonly,levels);
[NGTDMtextures] = getNGTDMtextures(NGTDM,countValid);
        t1 = string(struct2cell(textures));
        t2 = string(struct2cell(globaltexts));
        t3 = string(struct2cell(GLRLMtextures));
        t4 = string(struct2cell(GLSZMtextures));
        t5 = string(struct2cell(NGTDMtextures));
        towrite = [t1.',t2.',t3.',t4.',t5.'];
        allp{1,fxNum} = [{string(img.patientName.FamilyName)},{string(fractions(fxNum).name)},{'Equal'},{'64'},{string(tempstruct{1}.name)},towrite];
%         for r = 1:size(allp,2)
%          print{r} = char(allp(r));
%         end
end
        fx1print = allp{1,1};
        fx2print = allp{1,2};
        fx3print = allp{1,3};
        fx4print = allp{1,4};
        fx5print = allp{1,5};
        tprint = horzcat(fx1print,fx2print,fx3print,fx4print,fx5print);  
        for rr=1:length(tprint)
           print(1,rr) = {char(tprint(1,rr))}; 
        end
        xlsappend(xlfilename,print,1);
        system('taskkill /F /IM EXCEL.EXE');
  toc
% end
  %%%%%%%%%%%%%%%
