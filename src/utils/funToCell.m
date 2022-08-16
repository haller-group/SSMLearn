function dataOut = funToCell(dataIn,fun)
dataOut = dataIn;
for iTraj = 1:size(dataOut,1)
    dataOut{iTraj,2} = fun(dataIn{iTraj,2});
end
end