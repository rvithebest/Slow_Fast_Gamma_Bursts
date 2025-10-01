function [row,column,electrodeArray] = electrodePositionOnGrid(electrodeNum,gridType,subjectName,gridLayout)

if ~exist('subjectName','var');         subjectName=[];                 end
if ~exist('gridLayout','var');          gridLayout=2;                   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(gridType,'EEG')
    if gridLayout == 1 % EasyCap 64 electrodes
        
        electrodeArray = ...
            [00 00 00 00 01 61 02 00 00 00 00;
            00 00 53 00 39 00 40 00 54 00 00;
            00 11 47 03 33 17 34 04 48 12 00;
            29 55 25 41 21 20 22 42 26 56 30;
            00 13 49 05 35 18 36 06 50 14 00;
            31 57 27 43 23 62 24 44 28 58 32;
            00 15 51 07 37 19 38 08 52 16 00;
            00 00 59 00 45 63 46 00 60 00 00;
            00 00 00 00 09 64 10 00 00 00 00];
        
        %          electrodeLabelArray = ...
        %             ['   ' '   ' '   ' '   ' 'Fp1' 'Fpz' 'Fp2' '   ' '   ' '   ' '   ';
        %              '   ' '   ' 'AF7' '   ' 'AF3' '   ' 'AF4' '   ' 'AF8' '   ' '   ';
        %              '   ' 'F7' 'F5' 'F3' 'F1' 'Fz' 'F2' 'F4' 'F6' 'F8' '   ';
        %              'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'FC6' 'FT8' 'FT10';
        %              '   ' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' '   ';
        %              'TP9' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10';
        %              '   ' 'P7' 'P5' 'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' '   ';
        %              '   ' '   ' 'PO7' '   ' 'PO3' 'POz' 'PO4' '   ' 'PO8' '   ' '   ';
        %              '   ' '   ' '   ' '   ' 'O1' 'Oz' 'O2' '   ' '   ' '   ' '   '];
        
    elseif gridLayout == 2 % actiCap 64 electrodes
        electrodeArray = ...
            [00 00 00 01 00 00 00 02 00 00 00;
            00 33 00 34 00 00 00 35 00 36 00;
            00 03 37 04 38 05 39 06 40 07 00;
            41 42 08 43 09 00 10 44 11 45 46;
            00 12 47 13 48 14 49 15 50 16 00;
            17 51 18 52 19 53 20 54 21 55 22;
            00 23 56 24 57 25 58 26 59 27 00;
            28 60 00 61 00 62 00 63 00 64 32;
            00 00 00 29 00 30 00 31 00 00 00];
        
    elseif gridLayout == 3 % 10-20 system
        electrodeArray = ...
            [00 00 00 01 00 02 00 00 00 00;
            00 00 00 00 00 00 00 00 00 00;
            03 00 15 00 17 00 16 00 04 00;
            00 00 00 00 00 00 00 00 00 00;
            05 00 13 00 18 00 14 00 06 00;
            00 00 00 00 00 00 00 00 00 00;
            07 00 11 00 19 00 12 00 08 00;
            00 00 21 00 00 00 20 00 00 00;
            00 00 00 09 00 10 00 00 00 00];
        
        %          electrodeLabelArray = ...
        %             ['   ' '   ' '   ' 'Fp1' '   ' 'Fp2' '   ' '   ' '   ' '   ';
        %              '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ';
        %              'F7' '   ' 'F3' '   ' 'Fz' '   ' 'F4' '   ' 'F8' '   ';
        %              '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ';
        %              'T7' '   ' 'C3' '   ' 'Cz' '   ' 'C4' '   ' 'T8' '   ';
        %              '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ' '   ';
        %              'P7' '   ' 'P3' '   ' 'Pz' '   ' 'P4' '   ' 'P8' '   ';
        %              '   ' '   ' 'PO1' '   ' '   ' '   ' 'PO2' '   ' '   ' '   ';
        %              '   ' '   ' '   ' 'O1' '   ' 'O2' '   ' '   ' '   ' '   '];
    
    elseif gridLayout == 4  % coco EEG
        electrodeArray = ...
            [00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             10 00 00 00 00 00 00 00 11 ;
             00 07 00 06 00 08 00 09 00 ;
             03 00 02 00 01 00 04 00 05];
         
%          electrodeArray = ...           %20 Channel montage Alpa  
%             [00 17 00 00 00 00 00 18 00 ;
%              00 15 00 00 00 00 00 16 00 ;
%              00 00 00 00 00 00 00 00 00 ;
%              00 13 00 11 00 12 00 14 00 ;
%              00 00 00 00 02 00 00 00 00 ;
%              00 09 00 07 00 08 00 10 00 ;
%              05 00 03 00 01 00 04 00 06 ];
    
    elseif gridLayout == 5
        electrodeArray = ...
            [00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             00 00 00 00 00 00 00 00 00 ;
             106 00 00 00 00 00 00 00 107 ;
             00 103 00 102 00 104 00 105 00 ;
             99 00 98 00 97 00 100 00 101];
    end
    
    if electrodeNum<1 || electrodeNum>96
        disp('Electrode Number out of range');
        row=1;column=1;
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Microelectrode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(gridType,'Microelectrode')
    if strcmp(subjectName,'abu') || strcmp(subjectName,'rafiki')
        
        electrodeArray = ...
            [00 02 01 03 04 06 08 10 14 00;
            65 66 33 34 07 09 11 12 16 18;
            67 68 35 36 05 17 13 23 20 22;
            69 70 37 38 48 15 19 25 27 24;
            71 72 39 40 42 50 54 21 29 26;
            73 74 41 43 44 46 52 62 31 28;
            75 76 45 47 51 56 58 60 64 30;
            77 78 82 49 53 55 57 59 61 32;
            79 80 84 86 87 89 91 94 63 95;
            00 81 83 85 88 90 92 93 96 00];
        
    elseif strcmp(subjectName,'alpa') || strcmp(subjectName,'kesari') || strcmp(subjectName,'coco') || isempty(subjectName)
        
        electrodeArray = ...
            [00 88 78 68 58 48 38 28 18 00;
            96 87 77 67 57 47 37 27 17 08;
            95 86 76 66 56 46 36 26 16 07;
            94 85 75 65 55 45 35 25 15 06;
            93 84 74 64 54 44 34 24 14 05;
            92 83 73 63 53 43 33 23 13 04;
            91 82 72 62 52 42 32 22 12 03;
            90 81 71 61 51 41 31 21 11 02;
            89 80 70 60 50 40 30 20 10 01;
            00 79 69 59 49 39 29 19 09 00];
        
    elseif  strcmp(subjectName,'dona') || strcmp(subjectName,'hulk') || strcmp(subjectName,'jojo')
        
        electrodeArray = ...
            [96 90 84 78 72 66 60 54;
            95 89 83 77 71 65 59 53;
            94 88 82 76 70 64 58 52;
            93 87 81 75 69 63 57 51;
            92 86 80 74 68 62 56 50;
            91 85 79 73 67 61 55 49;
            00 00 00 00 00 00 00 00;
            48 42 36 30 24 18 12 06;
            47 41 35 29 23 17 11 05;
            46 40 34 28 22 16 10 04;
            45 39 33 27 21 15 09 03;
            44 38 32 26 20 14 08 02;
            43 37 31 25 19 13 07 01];
        
    elseif strcmp(subjectName,'tutu')|| strcmp(subjectName,'alpaH') || strcmp(subjectName,'kesariH')
        
        electrodeArray = ...
            [81 72 63 54 45 36 27 18 09;
            80 71 62 53 44 35 26 17 08;
            79 70 61 52 43 34 25 16 07;
            78 69 60 51 42 33 24 15 06;
            77 68 59 50 41 32 23 14 05;
            76 67 58 49 40 31 22 13 04;
            75 66 57 48 39 30 21 12 03;
            74 65 56 47 38 29 20 11 02;
            73 64 55 46 37 28 19 10 01;
            82 83 84 85 86 87 88 89 90];
    end
    
    if electrodeNum<1 || electrodeNum>96
        disp('Electrode Number out of range');
        row=1;column=1;
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ECoG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(gridType,'ECoG')
    
    if gridLayout==14
        goodElectrodeList = 1:4;
        electrodeArray = [1 2 3 4];
    elseif gridLayout==24
        goodElectrodeList = 1:8;
        electrodeArray = ...
            [1 2 3 4;
             5 6 7 8];
    elseif gridLayout==44
        goodElectrodeList = 1:16;
        electrodeArray = ...
            [01 02 03 04;
             05 06 07 08
             09 10 11 12
             13 14 15 16];
    else
        
        goodElectrodeList = [1:34 65:96];
        electrodeArray = ...
            [00 89 81 73 65 25 17 09 01 00;
            00 90 82 74 66 26 18 10 02 00;
            00 91 83 75 67 27 19 11 03 00;
            34 92 84 76 68 28 20 12 04 33;
            00 93 85 77 69 29 21 13 05 00;
            00 94 86 78 70 30 22 14 06 00;
            00 95 87 79 71 31 23 15 07 00;
            00 96 88 80 72 32 24 16 08 00];
    end
    if isempty(setdiff(goodElectrodeList,electrodeNum))
        disp('Electrode Number out of range');
    else
        [row,column] = find(electrodeArray==electrodeNum);
    end
end
end