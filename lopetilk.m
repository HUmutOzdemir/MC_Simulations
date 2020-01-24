%Number of molecules.
num_of_mols = 10000;
%Distance between release point and reciever(Micrometers).
distance = 16;
%Number of rows in the receiver group
num_of_rows = 150;
%Number of columns in the receiver group
num_of_columns = 15;
%I assume cells are square shaped. This is the length of an edge(Micrometers).
length_of_edge = 5;
%Diffusion coefficient(Micrometers^2/seconds).
coef = 79.4;
%Step size. (Seconds).
time_step = 10^-4;
%Total time. (Seconds).
time = 2;
%In micrometers per second
flow = [0, 0];
%Number of steps for a singal signal.
number_of_steps = 100;
%Molecule x, y, z coordinate matrix
mol_matrix = zeros(num_of_mols, 2);
%Standart deviation of the gaussian diffusion motion.
std = sqrt(2 * coef * time_step);
%Signal time of cells
signal_time_matrix = zeros(num_of_rows,num_of_columns);
%Signal step of cells
signal_step_matrix = zeros(num_of_rows,num_of_columns);

%List for bfs of signal.
bfs_cells = {};

%Simulation Loop
for i = 1:time/time_step
  
    %Generates  movements of current time step.
    mov_vectors = normrnd(0,std, size(mol_matrix));
    
    %Takes a copy of molecule matrix before calculations. Because if a
    %particle bumps it bounce back to its previous place. 
    prev_copy_mol_matrix = mol_matrix;
    
    %Calculates the new position of molecules by adding random movements
    %and current flow.
    mol_matrix = mol_matrix + mov_vectors + (flow * time_step);
    
    %Finds particles exceeds the channel
    exceeding_particles = find((mol_matrix(:,2)>=(num_of_rows*length_of_edge)/2) | (mol_matrix(:,2)<=((-1)*(num_of_rows*length_of_edge)/2)));
    
    %Bounce back to their previous position.
    mol_matrix(exceeding_particles,:) = prev_copy_mol_matrix(exceeding_particles,:);
    
    %Finds the molecules hitted to the wall of cells.
    hits = find(mol_matrix(:,1)>=distance);
    
    %Finds the molecules didn't hit to the wall of cells.
    not_hit = find(mol_matrix(:,1)<distance);
    
    %Finds the first hit cells.
    first_hit = mol_matrix(hits,:);
    first_hit= first_hit(:,2) + ((num_of_rows*length_of_edge)/2);
    first_hit = fix(first_hit/length_of_edge)+1;
    list = find(signal_time_matrix(first_hit,1)==0);
    
    %Saves the hitting step and time for first hit cells.
    signal_time_matrix(first_hit(list),1) = i*time_step;
    signal_step_matrix(first_hit(list),1) = number_of_steps; 
    %Adds the cell to the bfs list.
    temp_list = first_hit(list);
    for j=1:length(temp_list)
        bfs_cells = [bfs_cells,[temp_list(j),1]];
    end
    %Removes the hitted molecules
    mol_matrix = mol_matrix(not_hit,:);

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

freq = zeros(number_of_steps,1);
%Drop zeros in arrival times array.
arrival_times = signal_step_matrix(:);
res = arrival_times(arrival_times ~= 0);
%Creates a frequency array of the arriving steps.
for i = 1:length(res)
    freq(number_of_steps-res(i)+1) = freq(number_of_steps-res(i)+1) +1;
end
%Create a bar graph
figure(2);
bar(freq);