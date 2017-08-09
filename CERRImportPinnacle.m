% CERRImportPinnacle
% Import Native Pinnacle data into CERR plan format.
%
%APA, 08/28/2009
%
%
% Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team.
% 
% This file is part of The Computational Environment for Radiotherapy Research (CERR).
% 
% CERR development has been led by:  Aditya Apte, Divya Khullar, James Alaly, and Joseph O. Deasy.
% 
% CERR has been financially supported by the US National Institutes of Health under multiple grants.
% 
% CERR is distributed under the terms of the Lesser GNU Public License. 
% 
%     This version of CERR is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
% CERR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CERR.  If not, see <http://www.gnu.org/licenses/>.

%% --------------------------- Set Environment Vars for DCM4CHE
clear
flag = init_ML_DICOM;

if ~flag
    return;
end

%% --------------------------- Get the path of patient directory to be selected for import.
%pin_path = uigetdir(pwd','Select the directory containing native Pinnacle data');
pin_path = 'D:\User\Documents\MATLAB\Pinnacledata\ESO\Patient_9372';																 
% pin_path = 'D:\User\Documents\MATLAB\Pinnacledata\Patient_Lung';
if isnumeric(pin_path) || ~exist(fullfile(pin_path,'Patient'),'file')
    disp('Import of Pinnacle data aborted');
    return
end

tic;

CERRStatusString('Scanning native Pinnacle files');
%%--- get LA machine physic data
% phy_data = 'D:\用户目录\Documents\GitHub\Physics\ReadOnlyMachineDB';
% MachineList = {};
% Dirlist = dir(phy_data)
% linenum = 1;
% for line = 1:length(dir(phy_data))
%     if exist(fullfile(phy_data,Dirlist(line,1).name,'Pinnacle.Machine'),'file')
%         data_phys = read_1_pinnacle_file(fullfile(phy_data,Dirlist(line,1).name,'Pinnacle.Machine'),0);
%         MachineList{1,linenum} = data_phys;
%         linenum = linenum + 1;
%     end
% end
% machine_phy_data = 'D:\用户目录\Documents\GitHub\Pinnacle3_Native_parse\MachinePhyData.mat';
% machine_list = load(machine_phy_data,'-mat');
%machine_list.MachineList{1,N}.Name Varian 600CD
%machine_list.MachineList{1,N}.PhotonEnergyList.MachineEnergy.PhysicsData.O
%utpuFactor.SourceToCalibrationPointDistance 110
%machine_list.MachineList{1,N}.PhotonEnergyList.MachineEnergy.PhysicsData.O
%utpuFactor.DosePerMuAtCalibration 0.665
%% --------------------------- Read Patient data
data_pat = read_1_pinnacle_file(fullfile(pin_path,'Patient'),0);

%% ---------------------------  Populate indexS.scan field
%Use DCM4CHE to read scans
planC_w_scan = {};
for iScan = 1:length(data_pat.ImageSetList.ImageSet)
    try
        if length(data_pat.ImageSetList.ImageSet) == 1
            ImageName = data_pat.ImageSetList.ImageSet.ImageName;
        else
            ImageName = data_pat.ImageSetList.ImageSet{iScan}.ImageName;
        end
        scan_dcm_dir = [ImageName,'.DICOM'];
        planC_w_scan_tmp = readScanFromDir(fullfile(pin_path,scan_dcm_dir));
        %planC_w_scan_tmp =
        %readScanFromFile(fullfile(pin_path,scan_dcm_dir));

        %Populate fields like headInOut etc from Pinnacle file/s
        data_scan = read_1_pinnacle_file(fullfile(pin_path,[ImageName,'.header']),0);
        indexS_tmp = planC_w_scan_tmp{end};
        planC_w_scan_tmp{indexS_tmp.scan}.scanType = data_scan.modality;
        for iSInfo = 1:length(planC_w_scan_tmp{indexS_tmp.scan}.scanInfo)
            if data_scan.Version >= 9.0 %Pinnacle 9.0 dicom_image dicom begin using _x_start_dicom
                xOffset = data_scan.x_start_dicom + planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).grid2Units * (planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).sizeOfDimension2-1)/2;
                yOffset = data_scan.y_start_dicom + planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).grid1Units * (planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).sizeOfDimension1-1)/2;
            else %Pinnacle 8.X use x_start
                xOffset = data_scan.x_start + planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).grid2Units * (planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).sizeOfDimension2-1)/2;
                yOffset = data_scan.y_start + planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).grid1Units * (planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).sizeOfDimension1-1)/2;
            end
            planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).xOffset = xOffset;
            planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).yOffset = yOffset;
            planC_w_scan_tmp{indexS_tmp.scan}.scanInfo(iSInfo).headInOut = data_scan.patient_position;
        end
        if isempty(planC_w_scan)
            scanNum = 1;
        else
            scanNum = length(planC_w_scan) + 1;
        end

        planC_w_scan{scanNum} = planC_w_scan_tmp;

        db_nameC{scanNum} = data_scan.db_name;

    catch


    end
end

planC = initializeCERR;
indexS = planC{end};

for iScan = 1:length(planC_w_scan)
    planC_tmp = planC_w_scan{iScan};
    indexS_tmp = planC_tmp{end};
    %     planC{indexS.scan}(iScan) = planC_tmp{indexS_tmp.scan};
    planC{indexS.scan}=dissimilarInsert(planC{indexS.scan},planC_tmp{indexS_tmp.scan},iScan,[]);
end

%% --------------------------- Populate indexS.structures field
%loop over each plan and associate with a scan.
%Associated scan ID
% data_pat.PlanList.Plan{1}.PrimaryCTImageSetID

% BEGIN - Comment for structures
if length(data_pat.ImageSetList.ImageSet) == 1
    ImageSetIDv(iScan) = data_pat.ImageSetList.ImageSet.ImageSetID;
else
    for iScan = 1:length(data_pat.ImageSetList.ImageSet)
        ImageSetIDv(iScan) = data_pat.ImageSetList.ImageSet{iScan}.ImageSetID;
    end
end

%
plan_info = {};
for iPlan = 1:length(data_pat.PlanList.Plan)
    if length(data_pat.PlanList.Plan) == 1
        imageSetID = data_pat.PlanList.Plan.PrimaryCTImageSetID;
        planID = data_pat.PlanList.Plan.PlanID;
        plan_name  = data_pat.PlanList.Plan.PlanName;
        plan_dosimetrist = data_pat.PlanList.Plan.Dosimetrist;
        plan_physicist = data_pat.PlanList.Plan.Physicist;
        plan_comment = data_pat.PlanList.Plan.Comment;
        plan_TPS_version = data_pat.PlanList.Plan.PinnacleVersionDescription;
        plan_lock = data_pat.PlanList.Plan.PlanIsLocked;
    else
        imageSetID = data_pat.PlanList.Plan{iPlan}.PrimaryCTImageSetID;
        planID = data_pat.PlanList.Plan{iPlan}.PlanID;
        plan_name  = data_pat.PlanList.Plan{iPlan}.PlanName;
        plan_dosimetrist = data_pat.PlanList.Plan{iPlan}.Dosimetrist;
        plan_physicist = data_pat.PlanList.Plan{iPlan}.Physicist;
        plan_comment = data_pat.PlanList.Plan{iPlan}.Comment;
        plan_TPS_version = data_pat.PlanList.Plan{iPlan}.PinnacleVersionDescription;
        plan_lock = data_pat.PlanList.Plan{iPlan}.PlanIsLocked;
    end
     
    if (strfind(upper(plan_name),'QA') > 0)
        %data_pat.PlanList.Plan{iPlan} = {};
        continue;
    end
    
    if strfind(upper(plan_name),'COPY') %searching duplicate copy plan of previous plan, skip importing
        duplicate_plan = 0;
        for  indexPlan = 1:length(plan_list)
            if isempty(plan_list{indexPlan})
                continue;
            else
                ind_plan = plan_list{indexPlan}.index;
                if (imageSetID == plan_list{ind_plan}.imageSetID) && strcmp(plan_comment,plan_list{ind_plan}.plan_comment)
                    duplicate_plan = 1;
                end
            end
        end
        if duplicate_plan
            %data_pat.PlanList.Plan{iPlan} = {};
            continue;
        end
    end
    
    % if plan is not locked (plan_lock = 0), the plan parameter and plan
    % not trustable
    if null(plan_lock)
        continue;
    end
    %the QA plans do not contain CT images in the local plan dir, not parsing
    if (imageSetID == -1) 
        continue;
    end
    assocScanNum = find(ImageSetIDv == imageSetID);
    if length(planC{indexS.scan}) < assocScanNum
        continue;
    end   
    
    if isdir(fullfile(pin_path,['Plan_',num2str(planID)]))
        plan_list{iPlan}.index = iPlan;
        plan_list{iPlan}.imageSetID = imageSetID;
        plan_list{iPlan}.planID = planID;
        plan_list{iPlan}.plan_name = plan_name;
        plan_list{iPlan}.plan_dosimetrist = plan_dosimetrist;
        plan_list{iPlan}.plan_physicist = plan_physicist;
        plan_list{iPlan}.plan_comment = plan_comment;
        plan_list{iPlan}.plan_TPS_version = plan_TPS_version;
        plan_list{iPlan}.plan_lock = plan_lock;    
    end
    
    zValuesV = [planC{indexS.scan}(assocScanNum).scanInfo(:).zValue];
    
    %read structures from ROI files
    roi_filename = fullfile(pin_path,['Plan_',num2str(planID)],'plan.roi');
    data_plan_roi = read_1_pinnacle_file(roi_filename,0);
    if isfield(data_plan_roi,'roi') && ~isempty(data_plan_roi)
        for iStruct = 1:length(data_plan_roi.roi)
            %remove some auxilary contours (auto-create by some scripts)
            if ~isempty(regexp(upper(data_plan_roi.roi{iStruct}.name),'^1FS*', 'once')) 
                continue;
            end
            if ~isempty(regexp(upper(data_plan_roi.roi{iStruct}.name),'^1PTV*', 'once'))
                continue;
            end
            if ~isempty(regexp(upper(data_plan_roi.roi{iStruct}.name),'^PATIENT*', 'once'))
                continue;
            end
            if ~isempty(regexp(upper(data_plan_roi.roi{iStruct}.name),'^FS*', 'once'))
                continue;
            end
            planC{indexS.structures}(end+1) = newCERRStructure(assocScanNum,planC);
            planC{indexS.structures}(end).structureName = num2str(data_plan_roi.roi{iStruct}.name);
            planC{indexS.structures}(end).studyNumberOfOrigin = planID; %add planID mark to ROIs
            planC{indexS.structures}(end).roiNumber  = data_plan_roi.roi{iStruct}.num_curve; %contour numbers of current Roi
%             planC{indexS.structures}(end).structureDescription = data_plan_roi.roi{iStruct}.volume; %contour volume of current ROI
            if isfield(data_plan_roi.roi{iStruct},'curve')
                for iSlice = 1:length(data_plan_roi.roi{iStruct}.curve)
                    if iscell(data_plan_roi.roi{iStruct}.curve)
                        zVal = data_plan_roi.roi{iStruct}.curve{iSlice}.poin(1,3);
                    else
                        zVal = data_plan_roi.roi{iStruct}.curve.poin(1,3);
                    end
                    sliceNum = find((zValuesV-zVal).^2 < 1e-3);
                    if ~isempty(sliceNum) && isempty(planC{indexS.structures}(end).contour(sliceNum).segments(1).points)
                        if iscell(data_plan_roi.roi{iStruct}.curve)
                            planC{indexS.structures}(end).contour(sliceNum).segments(1).points = data_plan_roi.roi{iStruct}.curve{iSlice}.poin;
                        else
                            planC{indexS.structures}(end).contour(sliceNum).segments(1).points = data_plan_roi.roi{iStruct}.curve.poin;
                        end
                    elseif ~isempty(sliceNum)
                        if iscell(data_plan_roi.roi{iStruct}.curve)
                            planC{indexS.structures}(end).contour(sliceNum).segments(end+1).points = data_plan_roi.roi{iStruct}.curve{iSlice}.poin;
                        else
                            planC{indexS.structures}(end).contour(sliceNum).segments(end+1).points = data_plan_roi.roi{iStruct}.curve.poin;
                        end
                    end
                    switch upper(data_scan.patient_position)
                        case 'HFS' %+x,-y,-z
                            %planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,2) = 2*planC{indexS.scan}(assocScanNum).scanInfo(1).yOffset - planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,2);
                            %data(:,2) = 2*yOffset*10 - data(:,2);
                        case 'HFP' %-x,+y,-z
                            planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,1) = 2*planC{indexS.scan}(assocScanNum).scanInfo(1).xOffset - planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,1);
                            %data(:,1) = 2*xOffset*10 - data(:,1);
                        case 'FFS' %+x,-y,-z
                            %                             planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,2) = 2*planC{indexS.scan}(assocScanNum).scanInfo(1).yOffset - planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,2);
                            %data(:,2) = 2*yOffset*10 - data(:,2);
                        case 'FFP' %-x,+y,-z
                            planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,1) = 2*planC{indexS.scan}(assocScanNum).scanInfo(1).xOffset - planC{indexS.structures}(end).contour(sliceNum).segments(end).points(:,1);
                            %data(:,1) = 2*xOffset*10 - data(:,1);
                    end
                end
            end
        end
    end

    %Read transformation matrices
    volumeInfo_filename = fullfile(pin_path,['Plan_',num2str(planID)],'plan.VolumeInfo');
    if exist(volumeInfo_filename,'file')
        try
            data_plan_volumeInfo = read_1_pinnacle_file(volumeInfo_filename,0);
        catch
        end
    end

    %Scan Names
    for volNum = 1:length(data_plan_volumeInfo.VolumeDisplay)
        volNameC{volNum} = data_plan_volumeInfo.VolumeDisplay{volNum}.VolumeName;
    end
    for scanNum = 1:length(db_nameC)
        volNum = strcmpi(db_nameC{scanNum},volNameC);
        if ~isempty(volNum) % && data_plan_volumeInfo.VolumeDisplay{volNum}.Moveable
            planC{indexS.scan}(scanNum).transM = data_plan_volumeInfo.VolumeDisplay{volNum}.LocalToWorldTransform.Data.Points;
        end
    end   
   
end



% CERR Options
planC{indexS.CERROptions} = CERROptions;
planC = getRasterSegs(planC);
planC = setUniformizedData(planC);

% END - Comment for structures


%% ---------------------------  Populate indexS.dose field

for iPlan = 1:length(data_pat.PlanList.Plan)
   
    if length(data_pat.PlanList.Plan) == 1
        imageSetID = data_pat.PlanList.Plan.PrimaryCTImageSetID;
    else
        imageSetID = data_pat.PlanList.Plan{iPlan}.PrimaryCTImageSetID;
    end
    if (imageSetID == -1)
        continue;
    end
    assocScanNum = find(ImageSetIDv == imageSetID);

    if length(data_pat.PlanList.Plan) == 1
        planID = data_pat.PlanList.Plan.PlanID;
    else
        planID = data_pat.PlanList.Plan{iPlan}.PlanID;
    end
    fprintf('Plan_name: plan_%d\n',planID);
    planPath = fullfile(pin_path,['Plan_',num2str(planID)]);

    if exist(fullfile(planPath,'plan.Trial'),'file')
        plan_trial_data = read_1_pinnacle_file(fullfile(planPath,'plan.Trial'),0);
    else
        continue
    end
    
     % Read plan's points coordinate
    points_filename = fullfile(pin_path,['Plan_',num2str(planID)],'plan.Points');
    if exist(points_filename,'file')
        try
            data_plan_points = read_1_pinnacle_file(points_filename,0);
        catch
        end
    end
    
    % Read plan's machine infomation( physic data )
    machine_filename = fullfile(pin_path,['Plan_',num2str(planID)],'plan.Pinnacle.Machines');
    if exist(machine_filename,'file')
        try
            data_plan_machine = read_1_pinnacle_file(machine_filename,0);
        catch
        end
    end

    for iTrial = 1:length(plan_trial_data.Trial)

        if length(plan_trial_data.Trial) == 1
            Trial = plan_trial_data.Trial;
        else
            Trial = plan_trial_data.Trial{iTrial};
        end

        trialName = Trial.Name;
        fprintf('Trial_name: %s\n',trialName);

        plan = Trial;
        %         plan = read_pinnacle_plan_trial(fullfile(pin_path,['Plan_',num2str(planID)],'plan.Trial'),trialName);

        if isempty(plan)
            continue;
        end

        dosedim = [plan.DoseGrid__Dimension__Y plan.DoseGrid__Dimension__X plan.DoseGrid__Dimension__Z];

        if isempty(plan.BeamList)
            continue;
        end

        beams = plan.BeamList.Beam;

        prescriptions = plan.PrescriptionList.Prescription;
        if length(prescriptions) > 1
            N = length(prescriptions);
            for k = 1:N
                prenames{k} = prescriptions{k}.Name;
            end

            [selection,ok] = listdlg('ListString',prenames,'SelectionMode','single','Name','Prescription Selection','PromptString','Please select a prescription');
            if ok == 0
                return;
            end
            prescription = prescriptions{selection};
        else
            prescription = prescriptions;
            % 	selection = 1;
        end
        
        fprintf('Prescription = %s\n',prescription.Name);

        numbeam = length(beams);
        dosefilenames = find_pinnacle_dose_filenames(dosedim,planPath);
        if isempty(dosefilenames)
            continue;
        end
        dose = zeros(dosedim);
        for k = 1:numbeam
            beam = beams{k};
            fileName = beam.DoseVolume;
            indColon = findstr(fileName,':');
            fileName = fileName(indColon+1:end-1);
            fileName = sprintf('%3g',str2num(fileName));
            indZero = fileName == ' ';
            fileName(indZero) = '0';
            fileName = ['plan.Trial.binary.',fileName];
            if strcmpi(prescription.Name,beam.PrescriptionName) ~= 1
                continue;
            end
            fname = dosefilenames{k};
            %dose_per_beam = read_pinnacle_dose_file(fname,dosedim);
            

            try
                dose_per_beam = read_pinnacle_dose_file(fileName,dosedim,planPath);
                factor = beam.MonitorUnitInfo.NormalizedDose*beam.MonitorUnitInfo.CollimatorOutputFactor*beam.MonitorUnitInfo.TotalTransmissionFraction;
                machine_energy_list = data_plan_machine.Item_0.PhotonEnergyList.MachineEnergy;
                if length(machine_energy_list) == 1
                    dose_per_mu = machine_energy_list.PhysicsData.OutputFactor.DosePerMuAtCalibration;
                else
                    for line = 1:length(machine_energy_list)
                        if strcmp(machine_energy_list{line,1}.Name,beam.MachineEnergyName)
                            dose_per_mu = machine_energy_list{line,1}.PhysicsData.OutputFactor.DosePerMuAtCalibration;
                        end
                    end
                end                
                MUs = round(beam.MonitorUnitInfo.PrescriptionDose/(factor*dose_per_mu));
                
			catch
                %continue
            end

            fprintf('Beam %s, MUs = %d\n',beam.Name,MUs);
            %dose = dose + doses(:,:,:,k)*beam.ODM.MUForODM;
            %dose = dose + dose_per_beam*beam.ODM.MUForODM;
            dose = dose + dose_per_beam * beam.Weight;
        end

        %Normalize based on prescription
        %Structure Index associated with assocScanNum
        structNumV = getStructureAssociatedScan(1:length(planC{indexS.structures}),planC);
        structNumV = find(structNumV == assocScanNum);
        if isfield(prescription,'PrescriptionRoi')
            structIndex = strmatch(prescription.PrescriptionRoi,{planC{indexS.structures}(structNumV).structureName},'exact');
        elseif isfield(prescription,'PrescriptionPoint')
            structIndex = strmatch(prescription.PrescriptionPoint,{planC{indexS.structures}(structNumV).structureName},'exact');
        end
        if ~isempty(structIndex)
            structIndex = structIndex(1);
        end
        %Prescription Dose per fraction
        PrescriptionDose = prescription.PrescriptionDose;
        %Number of fractions
        numFractions = prescription.NumberOfFractions;
        %PrescriptionPercent 
        PrePercent = prescription.PrescriptionPercent * 0.01;
        %Total Dose
        totalPrescriptionDose = PrescriptionDose * numFractions / PrePercent; %Gy(perpercent*100)
        %Normalization method
        [jnk,NormalizationMethod] = strtok(prescription.NormalizationMethod);
        NormalizationMethod = deblank2(NormalizationMethod);

        %Populate planC dose field
        planC{indexS.dose}(end+1).doseArray = dose;
        clear dose;
        planC{indexS.dose}(end).doseType ='PHYSICAL';
        planC{indexS.dose}(end).doseUnits ='GRAYS';
        planC{indexS.dose}(end).doseScale =1;
        planC{indexS.dose}(end).fractionGroupID = trialName;
        planC{indexS.dose}(end).orientationOfDose ='TRANSVERSE';
        planC{indexS.dose}(end).numberOfDimensions =length(size(planC{indexS.dose}(end).doseArray));
        planC{indexS.dose}(end).sizeOfDimension1 = Trial.DoseGrid__Dimension__X;
        planC{indexS.dose}(end).sizeOfDimension2 = Trial.DoseGrid__Dimension__Y;
        planC{indexS.dose}(end).sizeOfDimension3 = Trial.DoseGrid__VoxelSize__Z;

        planC{indexS.dose}(end).horizontalGridInterval = Trial.DoseGrid__VoxelSize__X;
        planC{indexS.dose}(end).verticalGridInterval = - Trial.DoseGrid__VoxelSize__Y;
        planC{indexS.dose}(end).zValues = Trial.DoseGrid__Origin__Z:Trial.DoseGrid__VoxelSize__Z:(Trial.DoseGrid__Origin__Z + Trial.DoseGrid__VoxelSize__Z * (Trial.DoseGrid__Dimension__Z-1));
        planC{indexS.dose}(end).coord1OFFirstPoint = Trial.DoseGrid__Origin__X;
        planC{indexS.dose}(end).coord2OFFirstPoint = Trial.DoseGrid__Dimension__Y*Trial.DoseGrid__VoxelSize__Y + Trial.DoseGrid__Origin__Y;
        planC{indexS.dose}(end).doseUID = createUID('dose');
        try,  planC{indexS.dose}(end).assocScanUID = planC{indexS.scan}(assocScanNum).scanUID; end

        %Compute Normalization factor
        if ~isempty(structIndex)
            meanROIDose = meanDose(planC, structIndex, length(planC{indexS.dose}),'Absolute');
        end

        %Normalize Dose
        if ~isempty(structIndex)
            planC{indexS.dose}(end).doseArray = planC{indexS.dose}(end).doseArray * totalPrescriptionDose / meanROIDose;
        else
            planC{indexS.dose}(end).doseArray = planC{indexS.dose}(end).doseArray;
        end
        
        %calculate parameters
        for structNum = 1: length(planC{indexS.structures}) 
            if isempty(planC{indexS.structures}(1,structNum).roiNumber) %empty ROIs (with zero curves)
                continue;
            end
            if (planC{indexS.structures}(1,structNum).studyNumberOfOrigin ~= planID) %mul-plans, skip contours which not belong to this plan
                continue;
            end
            doseBinsV = [];
            volsHistV = [];
            [planC,doseBinsV,volsHistV] = getDVHMatrix(planC,structNum,length(planC{indexS.dose}));
            structVolume = calc_abs_vol(doseBinsV,volsHistV);
            maxdose = round(calc_maxDose(doseBinsV,volsHistV,2))*0.01;
            meandose = round(calc_meanDose(doseBinsV,volsHistV,2))*0.01;
            mindose = round(calc_minDose(doseBinsV,volsHistV,2))*0.01;
            fprintf('%d:%s,volume=%6.2f,maxdose=%6.2f,Meandose=%6.2f,mindose=%6.2f\n',structNum,planC{indexS.structures}(structNum).structureName,structVolume,maxdose,meandose,mindose);
            if strcmpi(planC{indexS.structures}(1,structNum).structureName,'LUNG_TOTAL')
                V5 = calc_Vx(doseBinsV, volsHistV, 500, 2)/structVolume;
                V20  = calc_Vx(doseBinsV, volsHistV, 2000, 2)/structVolume;
                V30  = calc_Vx(doseBinsV, volsHistV, 3000, 2)/structVolume;
                fprintf('Lung_toatal:V5=%3.2f,V20=%3.2f,V30=%3.2f\n',V5,V20,V30);
            end
        end

    end %trial

end %plan

%% ---------------------------  Populate indexS.beams field

toc;

 
% %% ---------------------------- Print Out ROIs' Statistical Info 
% indexS = planC{end};

% %% --------------------------- Save planC
% %save_planC(planC,planC{indexS.CERROptions});
% 
% %%--------------------------------- export to dicom-rt
destDir = 'D:\User\Documents\MATLAB\bb';
export_planC_to_DICOM(planC,destDir);


