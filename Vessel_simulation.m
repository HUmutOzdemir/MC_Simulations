%Number of molecules.
num_of_mols = 10000;
%Distance between release point and reciever(Micrometers).
distance = 10;
%Radius of cylinder channel(Micrometers).
channel_radius = 5;
%Coordinates of the cicular reciever(Micrometers).
reciever_coordinates = [distance, 0, 0];
%Diffusion coefficient(Micrometers^2/seconds).
coef = 79.4;
%Step size. (Seconds).
step = 10^-4;
%Total time. (Seconds).
time = 2;
%Flow vector in micrometers per second
flow = [0, 0, 0];
%Molecule x, y, z coordinate matrix
mol_matrix = zeros(num_of_mols, 3); 
%Arrival times of the molecules.
arrival_times = zeros(num_of_mols,1);
%Standart deviation of the gaussian diffusion motion.
std = sqrt(2 * coef * step);


%Simulation Loop
for i = 1:time/step
    
    %Generates  movements of current time step.
    mov_vectors = normrnd(0,std, size(mol_matrix));
    
    %Takes a copy of molecule matrix before calculations. Because if a
    %particle bumps it bounce back to its previous place. 
    prev_copy_mol_matrix = mol_matrix;
    
    %Calculates the new position of molecules by adding random movements
    %and current flow.
    mol_matrix = mol_matrix + mov_vectors + (flow * step);
    
    %Takes the y and z cordinates of all molecules for calculating
    %molecules that bounces its previous position.
    mol_matrix_yz = mol_matrix(:,2:3);
    
    %Squares y and z cordinates for all particles.
    mol_matrix_yz = mol_matrix_yz.^2;
    
    %Calculates the sum of squares for all rows(particles).
    sum_vector = sum(mol_matrix_yz,2);
    
    %Finds the particles which exceeds the channel. If sum of squares of y
    %and z cordinates is bigger than channel radius, then it is a exceeding
    %particle.
    %Exceeding Particles is a logical vector that holds 0 for no-exceeds, 1 for exceeds.
    exceeding_particles = find((arrival_times==0) & (sum_vector>channel_radius^2));
    
    %Bounce back to their previous position.
    mol_matrix(exceeding_particles,:) = prev_copy_mol_matrix(exceeding_particles,:);
    
    %Finds the number of molecules absorbed by receiver.
    %If the x cordinate is bigger than or equal to the distance between
    %releasing point and receiver, then it is absorbed.
    %Absorbs is a logical vector that holds 0 for no-absorbs, 1 for absorbs.
    absorbs = find((arrival_times==0) &  (distance<=mol_matrix(:,1)));
    
    %Store the arrival times of particles(absorb vector has value 1 in that row).
    arrival_times(absorbs) = i*step;
 
end

%Plotting the results.
%Drop initial zeros in arrival times array.
res = arrival_times(arrival_times ~= 0);


%Create X axis values and intervals for plotting.
xStart = 0;
dx = 0.01;
N = 200;
x = xStart + (0:N-1)*dx;

%Take histogram of values.
h = hist(res,N);
%Plot values
plot(x,h);



