% .........................................................................
% ****************  Communication II_Fall 2021_Dr.Emadi  ******************
% ******************************  HW-1  ***********************************
% ********************  DanesH Abdollahi - 9723053  ***********************
% .........................................................................
clc ; clear ; close all ;

%% Part 1
mu = 0 ; 
var = 1 ;% Source Variance
sigma = sqrt(var) ;
Ts = 0.0001 ;
x_min = mu - 10*sigma ;
x_max = mu + 10*sigma ;
x = x_min : Ts : x_max ;
pdf = ( 1/(sigma .* sqrt(2*pi)) ) .* exp( -0.5*(x-mu).^2 ./ var ) ; % Defining The Gaussian PDF

info = [x_min , x_max , Ts ] ;
a = [-5,-4,-2,0,1,3,5] ; % Quantizer Boundries
N = length(a) + 1 ;

x_hat = find_x_hat_For_Nonuniform_Quantizer(a , N , info , pdf ,x ) ; % Find Centroids Of Region's
Var_Of_Quantization_Noise =...
    Calculating_Distortion_For_RandomVar_Source( x , pdf , x_hat , a , info ) ;...
    % It's The same as Distortion

SQNR = 10*log10( var ./ Var_Of_Quantization_Noise ) ; % Calculating SQNR

% Print Outputs
fprintf("Part 1\n"); 
fprintf('The x_hat Matrix (Centroids of Regions) is :\n') ;
for i =1 : N
    fprintf("x_hat(%d) =%f   ",i,x_hat(i)) ;
end
fprintf('\n\nVariance Of Quantization Noise = %f\n',Var_Of_Quantization_Noise);
fprintf('SQNR = %f dB\n',SQNR);
fprintf('\n\n*******************************\n') ;


%% Part 2
% Part a & b
mu = 0 ; 
P = [1,5,10] ;
N_main = [4,8,16,32] ;
fprintf("Part 2\n");
for z = 1 : length(P)
    var = P(z) ; %  Source Variance
    sigma = sqrt(var) ;
    Ts = 0.001 ;
    x_min = mu - 10*sigma ;
    x_max = mu + 10*sigma ;
    x = x_min : Ts : x_max ;
    pdf = ( 1/(sigma .* sqrt(2*pi)) ) .* exp( -0.5*(x-mu).^2 ./ var ) ; % Defining The Gaussian PDF
    info = [x_min , x_max , Ts ] ;
    epsilon = Ts ;
    SQNR_For_plot = zeros(1,length(N_main));
    for j = 1 : length(N_main)
        N = N_main(j) ;
        
        x_hat = sort( x_min + rand(1,N) .* (x_max - x_min) ) ; ...
            % Generating N Random Numbers in period (x_min : x_max) for Initial x_hat
        
        g=1 ;
        while 1       % Calculating a(The Boundries) & x_hat By Max_Lloyd Algotithm
            a = Find_a_For_NonUniform_Quantizer (x_hat , N) ;
            if ( g>1 ) && ( norm(a - a_pre) < epsilon ) % Termination condition of the loop
                Iteration_numb = g ;
                break ;
            end
            x_hat = find_x_hat_For_Nonuniform_Quantizer (a , N , info , pdf ,x ) ;
            a_pre = a ;
            g = g+1 ;
        end

        Distortion = Calculating_Distortion_For_RandomVar_Source( x , pdf , x_hat , a , info ) ;...
            % Calculating Distortion
        SQNR = 10*log10( var ./ Distortion ) ; % Calculating SQNR
        
        SQNR_For_plot(j) = SQNR ;
        % Outputs
        fprintf("N = %d\n",N) ; % Print The Ouput's
        fprintf("Var = %d\n",var) ;
        fprintf("Number Of Iterations = %d\n\n", Iteration_numb) ;
        
        for i =1 : N-1
            fprintf("a(%d) =%f   ",i,a(i)) ;
        end
        fprintf("\n\n");
        
        for i =1 : N
            fprintf("x_hat(%d) =%f   ",i,x_hat(i)) ;
        end
        fprintf("\n\n");
        
        fprintf("Distortion = %f\n", Distortion) ;
        
        fprintf("SQNR = %f dB\n", SQNR) ;
        fprintf("\n\n************************************************************************\n\n") ;
  
    end
    figure(z);
    stem(N_main , SQNR_For_plot,'*','LineWidth',2) ;
    title(sprintf('SQNR For a Gaussian Zero-mean R.V with Variance %d', P(z)));
    xlabel('Number of Quantization Levels');
    ylabel('SQNR (dB)') ;
    grid on ;
      
end


%% Functions

function y = integ(x , info , a , b) % The Func. For Calculating The integral Of x-function form a to b
    if a > 0
        a = fix((a+info(3)/2)/info(3)) * info(3) ;
    elseif a < 0
        a = fix((a-info(3)/2)/info(3)) * info(3) ;
    elseif a == 0 
        a = fix((a/2)/info(3)) * info(3) ;
    end
    
    if b > 0
        b = fix((b+info(3)/2)/info(3)) * info(3) ;
    elseif b < 0
        b = fix((b-info(3)/2)/info(3)) * info(3) ;
    elseif b == 0 
        b = fix((b/2)/info(3)) * info(3) ;
    end
    
    min = ((a - info(1))/info(3)) + 1 ;
    max = ((b - info(1))/info(3)) + 1 ;
    y = 0 ;
    for i = ceil(min) :  floor(max) 
        y = y + ( x(i) .* info(3) );
    end
    
end
% ------------------------------------------------------------------------------------------
function D = Calculating_Distortion_For_RandomVar_Source( x , pdf , x_hat , a , info )
    N = length(x_hat) ;
    D = 0 ; 
    D = D + integ( ((x - x_hat(1)).^2) .* pdf , info , info(1) , a(1) ) ;
    D = D + integ( ((x - x_hat(N)).^2) .* pdf , info , a(N-1) , info(2) ) ;
    
    for i=2 : N-1
        D = D + integ( ((x - x_hat(i)).^2) .* pdf  , info , a(i-1) , a(i) ) ;
    end


end
% ------------------------------------------------------------------------------------------------------------
function y = find_x_hat_For_Nonuniform_Quantizer(a , N , info , pdf , x) % The Func. For Calcultaing x_hat(i)

    x_hat(1) = integ(pdf .* x , info , info(1) , a(1) ) ./ integ(pdf  , info , info(1) , a(1) ) ;
    x_hat(N) = integ(pdf .* x , info , a(N-1) , info(2) ) ./ integ(pdf  , info , a(N-1) , info(2) ) ;

    for i=2 : N-1
        x_hat(i) = integ( pdf .* x , info , a(i-1) , a(i)) ./ integ(pdf  , info , a(i-1) , a(i) ) ;
    end
    y = x_hat ;
end
% -------------------------------------------------------------------------------------------------------------
function a = Find_a_For_NonUniform_Quantizer(x_hat , N) % The Func. For Calcultaing a(i)
    a = zeros(1,N-1);
    for i=1 : N-1
        a(i) = 0.5 .* ( x_hat(i) + x_hat(i+1) );
    end
end
