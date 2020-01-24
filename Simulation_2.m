%Number of molecules.
num_of_mols = 10000;
%Distance between release point and reciever(Micrometers).
distance = 48;
%Number of rows in the receiver group
num_of_rows = 1000;
%Number of columns in the receiver group
num_of_columns = 500;
%I assume cells are square shaped. This is the length of an edge(Micrometers).
radius_of_cell = 2.5;
%Diffusion coefficient(Micrometers^2/seconds).
coef = 79.4;
%Step size. (Seconds).
time_step = 10^-4;
%Total time. (Seconds).
time = 2;
%In micrometers per second
flow = [0, 0];
%Maximum number of Calcium Signalling Steps
max_steps = inf;
%Probability of Calcium Signaling
probability_of_signal = 0.7;
%Probability decrease in every Signal between cells
probability_decrease = 0.1;
%Probability increase in every time step.
probability_increase = 0.03;
%Molecule x, y, z coordinate matrix
mol_matrix = zeros(num_of_mols, 2);
%Standart deviation of the gaussian diffusion motion.
std = sqrt(2 * coef * time_step);
%Signal time of cells
signal_time_matrix = zeros(num_of_rows,num_of_columns);
%Signal step of cells
signal_step_matrix = zeros(num_of_rows,num_of_columns);
%Signalling probability of cells.
signal_probability_matrix = zeros(num_of_rows,num_of_columns);

%List for bfs of signal.
bfs_cells = {};

%Simulation Loop
for i = 1:time/time_step
    
     %To store new cells will make calcium signaling.
    temp = {};
    %Calcium Signaling Loop
    for j = 1:length(bfs_cells)
        count = 0;
        %Coordinates fo current cell in the cell matrix.
        current = cell2mat(bfs_cells(j));
        x = current(1);
        y = current(2);
        current_step = signal_step_matrix(x,y);
        %Max number of steps check
        if current_step == max_steps
            continue;
        end
        %Simulates the Calcium Signaling as BFS. A cell signals its
        %right, left, top and bottom neighboors if signal is not completed
        %and target cell was not signalled.
        if x+1<=num_of_rows && signal_step_matrix(x+1,y)==0 && rand <= signal_probability_matrix(x,y)
            signal_step_matrix(x+1,y) = current_step+1;
            signal_time_matrix(x+1,y) = i*time_step;
            signal_probability_matrix(x+1,y) = signal_probability_matrix(x,y) - probability_decrease;
            temp = [temp,[x+1,y]];
        else
            count = count+1;
        end
        if x-1>0 && signal_step_matrix(x-1,y)==0 && rand <= signal_probability_matrix(x,y)
            signal_step_matrix(x-1,y) = current_step+1;
            signal_time_matrix(x-1,y) = i*time_step;
            signal_probability_matrix(x-1,y) = signal_probability_matrix(x,y) - probability_decrease;
            temp = [temp,[x-1,y]];
        else
            count = count+1;
        end
        if y+1<=num_of_columns && signal_step_matrix(x,y+1)==0 && rand <= signal_probability_matrix(x,y)
            signal_step_matrix(x,y+1) = current_step+1;
            signal_time_matrix(x,y+1) = i*time_step;
            signal_probability_matrix(x,y+1) = signal_probability_matrix(x,y) - probability_decrease;
            temp = [temp,[x,y+1]];
        else
            count = count+1;
        end
        if y-1>0 && signal_step_matrix(x,y-1)==0 && rand <= signal_probability_matrix(x,y)
            signal_step_matrix(x,y-1) = current_step+1;
            signal_time_matrix(x,y-1) = i*time_step;
            signal_probability_matrix(x,y-1) = signal_probability_matrix(x,y) - probability_decrease;
            temp = [temp,[x,y-1]];  
        else
            count = count+1;
        end
        %if there is a cell that can be signalled
        if count ~= 4
            signal_probability_matrix(x,y) = signal_probability_matrix(x,y) + probability_increase;
            temp = [temp,[x,y]];
        end
    end
    %Updates calcium signaling iteration for next time step.
    bfs_cells = temp;

    %Generates  movements of current time step.
    mov_vectors = normrnd(0,std, size(mol_matrix));
    
    %Takes a copy of molecule matrix before calculations. Because if a
    %particle bumps it bounce back to its previous place. 
    prev_copy_mol_matrix = mol_matrix;
    
    %Calculates the new position of molecules by adding random movements
    %and current flow.
    mol_matrix = mol_matrix + mov_vectors + (time_step*flow);
    
    %Finds particles exceeds the channel
    exceeding_particles = find((mol_matrix(:,2)>=(num_of_rows*radius_of_cell) | (mol_matrix(:,2)<=((-1)*(num_of_rows*radius_of_cell)))));
    
    %Bounce back to their previous position.
    mol_matrix(exceeding_particles,:) = prev_copy_mol_matrix(exceeding_particles,:);
    
    %Finds the molecules hitted to the wall of cells.
    possible_hits = find(mol_matrix(:,1)>=distance);
    
    %Finds the first hit cells.
    possible_hits_pos = mol_matrix(possible_hits,:) + [0 (num_of_rows*radius_of_cell)];
    which_cell = possible_hits_pos(:,2);
    which_cell = fix(which_cell/(2*radius_of_cell))+1;
    cell_y = (2*which_cell - 1)*radius_of_cell;
    possible_hits_pos = possible_hits_pos - [(distance+radius_of_cell)*ones(length(cell_y),1) , cell_y];
    possible_hits_pos = sum(possible_hits_pos.^2,2);
    hits = find(possible_hits_pos<=radius_of_cell^2);
    list = find(signal_step_matrix(which_cell(hits),1)==0);
    
    temp = which_cell(hits);
    %Saves the hitting step and time for first hit cells.
    signal_time_matrix(temp(list),1) = i*time_step;
    signal_step_matrix(temp(list),1) = 1;
    signal_probability_matrix(temp(list),1) = probability_of_signal;
    %Adds the cell to the bfs list.
    for j=1:length(temp)
        bfs_cells = [bfs_cells,[temp(j),1]];
    end
    %Removes the molecules that were hit.
    mol_matrix(possible_hits(hits),:) = [];

end

%Drop zeros in arrival times array.
arrival_times = signal_time_matrix(:);
res = arrival_times(arrival_times ~= 0);
%Create X axis values and intervals for plotting.
xStart = 0;
dx = 0.01;
N = 200;
x = xStart + (0:N-1)*dx;
%Take histogram of values.
h = hist(res,N);
%Plot values
figure(1);
plot(x,h);


%Drop zeros in arrival times array.
arrival_steps = signal_step_matrix(:);
res2 = arrival_steps(arrival_steps ~= 0);
%Creates a frequency array of the arriving steps.
m = max(res2);
freq = zeros(m,1);
for i = 1:length(res2)
    freq(res2(i)) = freq(res2(i)) +1;
end
%Create a bar graph
figure(2);
bar(freq);