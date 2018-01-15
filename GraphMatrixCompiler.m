function [ Solution1, Solution2] = GraphMatrixCompiler( Matrix )
%Takes a Matrix for a graph, with the first column as the x axis and the rest of the columns
%as the other variable and turns the matrix into an array so that the
%Matlab plotter can plot it correctly. Solution1 is the main solution and
%solution 2 holds any seperate solution.

regions = 0; %counter or stability regions
PreNumSol = 0; %represents number of solutions in last row. Will be used later
Solution2 = [];
%Step 1, initial scan Matrix to count stability regions
for i=1:length(Matrix(:,1))
    NewNumSol=1; %counts the number of solutions in this region
    for j=2:length(Matrix(i,:)) %find number of solutions in this region
        if Matrix(i,j) ~= 0
            NewNumSol= NewNumSol + 1;
        end
    end
    if NewNumSol ~= PreNumSol %if this row has a different number of solutions than the last then it is a new region
        regions= regions+1;
    end
    PreNumSol=NewNumSol;
end

%Step 2, scan Matrix again and find start and end sites of each region
Tracker = zeros(regions, 3); %tracks the start and end point of each region and the number of solutions in the region
PreNumSol = 0;
Tracker(1,1)=1; %first region will always have a start point of 1
Tracker(regions,2)= length(Matrix(:,1)); %last region will always have the last end point
Counter = 1; %New counter to track region number
for i=1:length(Matrix(:,1))
    NewNumSol=0; %counts the number of solutions in this region
    for j=2:length(Matrix(i,:)) %find number of solutions in this region
        if Matrix(i,j) ~= 0
            NewNumSol= NewNumSol + 1;
        end
    end
    if NewNumSol ~= PreNumSol %if this row has a different number of solutions than the last then it is a new region
        if PreNumSol == 0
            Tracker(Counter,3)=NewNumSol;
        else
            Tracker(Counter,2) = i-1; %The last point was the end of the last region
            Counter = Counter+1;
            Tracker(Counter,1) = i; %This i value is the start of the new region
            Tracker(Counter,3) = NewNumSol; %number of solutions from new region
        end
     end
    PreNumSol=NewNumSol;
end

%Step 3, Create a Matrix representing seperate arrays of every solution
%within every section
TotSolutions = 0; %count pointer for total solutions segments
for i=1:length(Tracker(:,1)) %count total segments
    TotSolutions= TotSolutions+Tracker(i,3);
end

SolutionMatrix = zeros(TotSolutions, 5); % SolutionMatrix represents all solution segments
Counter=1;
for i=1:length(Tracker(:,1)) %Generate filled in Solution Matrix
    for j=1:Tracker(i,3)
        SolutionMatrix(Counter,1)=Tracker(i,1); %Start matrix i value
        SolutionMatrix(Counter,2)=Matrix(Tracker(i,1),j+1); %Start Y value
        SolutionMatrix(Counter,3)=Tracker(i,2); %End matrix i value
        SolutionMatrix(Counter,4)=Matrix(Tracker(i,2),j+1); %End Y value
        SolutionMatrix(Counter,5)=j+1; %Column within original Matrix
        Counter = Counter+1;
    end
end
% SolutionMatrix
% Step 4, Compile the sections for a solution
Solution1 = Matrix(1:SolutionMatrix(1,3),[1,2]);
SolutionMatrix(1,:) = []; %removes the segment from the solution matrix
StartEnd = -1; % Indicates if the  point is a start point or an end point. the direction of the algorithm. 1 = end,forwards, 0 = start, backwards, -1 = no point 
SegmentPointer = [1, 100]; % Indicates the best candidate for the next segment and the difference between the last point and the candidate
Counter=0; %counts the number of passes through without segment compilation. If it equals 2 then we must find a second solution
while isempty(SolutionMatrix)==false
    SegmentPointer = [];
    for i=1:length(SolutionMatrix(:,1))
        if (Matrix(SolutionMatrix(i,1),1) - Solution1(end,1) == 0) || abs(abs(Matrix(SolutionMatrix(i,1),1) - Solution1(end,1)) - abs(Solution1(length(Solution1)-1,1)-Solution1(end,1))) <= 2e-2 %if the start point of solution i is close in x value
            if isempty(SegmentPointer)
                SegmentPointer= [i,abs(SolutionMatrix(i,2) - Solution1(end,2))];
                StartEnd = 1;
            end
            if abs(SolutionMatrix(i,2) - Solution1(end,2)) <= SegmentPointer(2)
                SegmentPointer= [i,abs(SolutionMatrix(i,2) - Solution1(end,2))];
                StartEnd = 1;
            end
        end
        
        if Matrix(SolutionMatrix(i,3),1) - Solution1(end,1) == 0 || abs(abs(Matrix(SolutionMatrix(i,3),1) - Solution1(end,1)) - abs(Solution1(length(Solution1)-1,1)-Solution1(end,1))) <= 2e-2 %if the start point of solution i is close in x value
            if isempty(SegmentPointer)
                SegmentPointer= [i,abs(SolutionMatrix(i,4) - Solution1(end,2))];
                StartEnd = 0;
            end
            if abs(SolutionMatrix(i,4) - Solution1(end,2)) <= SegmentPointer(2)
                SegmentPointer= [i,abs(SolutionMatrix(i,4) - Solution1(end,2))];
                StartEnd = 0;
            end
        end
    end
    
    %Compile the best segment into the Solution
    if StartEnd == -1 %if no segment was found to attach to the end of Solution1, flip matrix to run in opposite direction
        Solution1 = flipud(Solution1);
        Counter = Counter+1;
    end
    if StartEnd == 1 %Compile into solution in the forward direction
        Solution1 = [Solution1; Matrix(SolutionMatrix(SegmentPointer(1),1):SolutionMatrix(SegmentPointer(1),3),[1,SolutionMatrix(SegmentPointer(1),5)])];
        SolutionMatrix(SegmentPointer(1),:) = []; %removes the segment from the solution matrix
        StartEnd = -1;
    end
    if StartEnd == 0 %Compile into solution in
        Solution1 = [Solution1; flipud(Matrix(SolutionMatrix(SegmentPointer(1),1):SolutionMatrix(SegmentPointer(1),3),[1,SolutionMatrix(SegmentPointer(1),5)]))];
        SolutionMatrix(SegmentPointer(1),:) = []; %removes the segment from the solution matrix
        StartEnd = -1;
    end
    if Counter ==2
        %Solution should only have to flip once so that it runs completely in
        %the forward and backwards direction. On the second count it is fliped
        %back to the origional direction. If there are any solutions left in
        %the Solution Matrix they must be discontinious from the main solution
        %This step generates a new matrix only containing the left over
        %solutions and recalls this function with the new matrix for Solution2
        Column= 0; % holds the number of columns for the new matrix
        for i=1:length(SolutionMatrix(:,1))
            if SolutionMatrix(i,5) > Column
                Column = SolutionMatrix(i,5);
            end
        end
        NewMatrix = zeros(SolutionMatrix(end,3)-SolutionMatrix(1,1)+1,Column); %New matrix for solution2
        NewMatrix(:,1) = Matrix(SolutionMatrix(1,1):SolutionMatrix(end,3),1); %Generate X values for New Matrix
        %Now fill in the Y values
        for i=1:length(SolutionMatrix(:,1))
            NewMatrix((SolutionMatrix(i,1)-SolutionMatrix(1,1)+1):(SolutionMatrix(i,3)-SolutionMatrix(1,1)+1),SolutionMatrix(i,5)) = ...
                Matrix(SolutionMatrix(i,1):SolutionMatrix(i,3),SolutionMatrix(i,5));
        end
        %New Code
        for i=length(NewMatrix(1,:)):-1:1 %New Matrix should not have any zero filled arrays
            if any(NewMatrix(:,i)) == 0
                NewMatrix(:,i)= [];
            end
        end
        % End of New Code
        Solution2 = GraphMatrixCompiler(NewMatrix);
        return
    end
end
if Solution1(1,1) == Solution1(end,1) %solution is most likely a circle and so the two lines will be merged
    Solution1 = [Solution1; Solution1(1,1), Solution1(1,2)];
end
end