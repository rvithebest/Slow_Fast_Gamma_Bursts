function distance = getElectrodeDistance(lfpElectrodes,spkElectrode)
gridType='Microelectrode';
[spkRow,spkColumn,electrodeArray] = electrodePositionOnGrid(spkElectrode,gridType,'alpaH');
distance = zeros(1,length(lfpElectrodes));
for i=1:length(lfpElectrodes)
    [lfpRow,lfpColumn] = find(electrodeArray==lfpElectrodes(i));
    distance(i) = 400*sqrt((lfpRow-spkRow).^2 + (lfpColumn-spkColumn).^2);
end
end