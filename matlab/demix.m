function demix(outPath, path, specStart, numJobs, numScans,algName, lambda1, lambda2, alpha, deisotope, calcPrecursorMass, globalTol, copy)
  % setup cvx
  install_cvx()

  filename = strcat(num2str(specStart), '.mgf');
  fileID = fopen(strcat(outPath, '/', filename), 'w');
  %Rout = fopen(strcat(outPath, '/', num2str(specStart), '.tab'), 'w');
  Rout = 1
  for i=specStart:numJobs:numScans
        if exist(strcat(path, num2str(i), '.tab'), 'file')
	  demix_spectrum(path, fileID, num2str(i), algName, lambda1, lambda2, alpha, deisotope, calcPrecursorMass, globalTol, copy, Rout);
        end
    end

  fclose(fileID);
  %fclose(Rout);
end

function demix_spectrum(path, outFile, expName, algName, lambda1, lambda2, alpha, deisotope, calcPrecursorMass, globalTol, copy, Rout)
    % read in vector b
    fileID = fopen(strcat(path, 'b_', expName, '.bin'));
    b = (fread(fileID, 'double'));
    fclose(fileID);
    minObs = min(b(b>0))-1;
    %b(b>0) = b(b>0) - minObs;
    n = size(b, 1);
    % read in matrix A
    fileID = fopen(strcat(path, 'A_', expName, '.bin'));
    A = (fread(fileID, 'double'));
    fclose(fileID);
    A = reshape(A, n, []);
    m = size(A,2);

    % read in start/end indices of precursor options in A
    fileID = fopen(strcat(path, 'indices_', expName, '.bin'));
    indices = fread(fileID, 'int32');
    indices = reshape(indices,2,[])';
    fclose(fileID);
    % read in group weights of precursor options in A
    %fileID = fopen(strcat(path, 'groupWeights_', expName, '.bin'));
    %groupWeights = 1./fread(fileID, 'double');
    %fclose(fileID);
    % read in precursor option titles
    precursorOptions = importdata(strcat(path, 'precursorOptions_', expName, '.tab'));
    % read in scan details
    scanDetails = importdata(strcat(path, expName, '.tab'));
    % read in mzValues
    mzValues = importdata(strcat(path, 'mz_', expName, '.tab'));

    if size(indices,1) > 1
    %if size(indices,1) <= 1
    %    [x, cvx_status] = NNLS(A, b, m);
        if strcmp(algName, 'OLS')
            [x, cvx_status] = OLS(A, b, m);
        elseif strcmp(algName, 'NNLS')
            [x, cvx_status] = NNLS(A, b, m, indices, precursorOptions, 2.0, lambda1);  
        elseif strcmp(algName, 'NNLS-A')
            [x, cvx_status] = NNLS_A(A, b, m, lambda1, precursorOptions, indices);
        elseif strcmp(algName, 'NNLS-L1')
            [x, cvx_status] = NNLS_L1(A, b, m, lambda1, indices, precursorOptions, 2.0);
        elseif strcmp(algName, 'NNLS-wL1')
            [x, cvx_status] = NNLS_weighted_L1(A, b, m, groupWeights, indices, lambda1);
        elseif strcmp(algName, 'NNLS-wL1-wL2')
            [x, cvx_status] = NNLS_weighted_L1_L2(A, b, m, groupWeights, indices, lambda1, lambda2);
        elseif strcmp(algName, 'NNLS-Linf')
            [x, cvx_status] = NNLS_group_inf(A, b, m, groupWeights, indices, lambda1);
        elseif strcmp(algName, 'NNLS-wLinf')
            [x, cvx_status] = NNLS_weight_group_inf(A, b, m, groupWeights, indices, lambda1);
        elseif strcmp(algName, 'NNLS-L1-Linf')
            [x, cvx_status] = NNLS_L1_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2);
        elseif strcmp(algName, 'NNLS-wL1-Linf')
            [x, cvx_status] = NNLS_weighted_L1_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2);
        elseif strcmp(algName, 'NNLS-L1-wLinf')
            [x, cvx_status] = NNLS_L1_weighted_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2);
        elseif strcmp(algName, 'NNLS-wL1-wLinf')
            [x, cvx_status] = NNLS_weighted_L1_weighted_group_inf(A, b, m, groupWeights, indices, lambda1);
        elseif strcmp(algName, 'NNLS-sparseGroupLasso')
            [x, cvx_status] = NNLS_sparse_group_lasso(A, b, m, n, indices, lambda1, lambda2, alpha, precursorOptions, 2);
        elseif strcmp(algName, 'NNLS-sparseGroupLasso2')
            [x, cvx_status] = NNLS_sparse_group_lasso2(A, b, m, n, indices, lambda1, lambda2, alpha, precursorOptions, 2);
        end
    end

    %if size(indices,1) == 1
    %    [x, cvx_status] = NNLS(A, b, m);
    %end

    if size(indices,1) > 1 && strcmp(cvx_status, 'Solved')
    %if strcmp(cvx_status, 'Solved') 
        if deisotope
            spectra = computeMonoSpectra(A, x, indices);
        else
            spectra = computeSpectra(A, x, indices);
        end

        writeMGF(spectra, mzValues, precursorOptions, calcPrecursorMass, scanDetails, outFile, globalTol, copy, b, minObs);

        %plotSpectra(spectra, b, precursorOptions, algName, expName, outPath);
        %writeR(spectra, mzValues, precursorOptions, scanDetails, Rout, b);

        %monoSpectra = computeMonoSpectra(A, x, indices);
        %spectra = computeSpectra(A, x, indices);
            
        %writeMGF2(spectra, monoSpectra, mzValues, precursorOptions, calcPrecursorMass, scanDetails, outFile, globalTol, copy, b, minObs); 


    else
        writeOriginalMGF(b, mzValues, scanDetails, outFile, globalTol);
    end

    %sprintf("Conditioning: ")
    %dA = decomposition(A);
    %isIllConditioned(dA)
    %rank(dA)
    %min(size(A))
    %cond(A)

end


function install_cvx()
    currentDir = pwd;
    cd /nas/longleaf/home/dennisg/cvx/cvx
    cvx_setup
    cd(currentDir)
end

function [x, cvx_status] = OLS(A, b, m)
    cvx_begin
        variable x(m);
        minimize( norm(A*x-b) );
    cvx_end
end


function [x, cvx_status] = NNLS(A, b, m, indices, precursorOptions, fudge, lambda)
    cvx_begin
        variable x(m);
        minimize( norm(A*x-b) );
        subject to
            x >= 0;
    cvx_end
    sprintf("TIC: %f", sum(A*x)/sum(b))
end


% NNLS + L1 penalty
function [x, cvx_status] = NNLS_L1(A, b, m, lambda, indices, precursorOptions, fudge)

    lambda = lambda / 2^(size(indices,1) - 2); 

    cvx_begin
        variable x(m);
        minimize( norm(A*x-b) + (lambda * norm(x,1)));
        subject to
            x >= 0; 
    cvx_end

    sprintf("TIC: %f", sum(A*x)/sum(b))
end


% NNLS weighted L1 penalty
function [x, cvx_status] = NNLS_weighted_L1(A, b, m, groupWeights, indices, lambda)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            i1 = indices(i,1)+1;
            i2 = indices(i,2)+1;
            xi = x(i1:i2);
           
            objective = objective + (lambda * groupWeights(i) * norm(xi, 1));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end

% NNLS weighted L1 penalty + l2 penalty
function [x, cvx_status] = NNLS_weighted_L1_L2(A, b, m, groupWeights, indices, lambda1, lambda2)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            objective = objective + (groupWeights(i) * lambda1 * norm(x(indices(i,1)+1:indices(i,2)+1), 1)) + (lambda2 * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), 2));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end

% NNLS group infinity norm
function [x, cvx_status] = NNLS_group_inf(A, b, m, groupWeights, indices, lambda)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            objective = objective + (lambda * norm(x(indices(i,1)+1:indices(i,2)+1), inf));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end


% NNLS weighted group infinity norm
function [x, cvx_status] = NNLS_weight_group_inf(A, b, m, groupWeights, indices, lambda)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            objective = objective + (lambda * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), inf));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end


% NNLS group infinity norm + L1 penalty
function [x, cvx_status] = NNLS_L1_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b) + (lambda1 * norm(x,1));
        for i=1:size(groupWeights,1)
            objective = objective + (lambda2 * norm(x(indices(i,1)+1:indices(i,2)+1), inf));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end


% NNLS group infinity norm + weighted L1 penalty
function [x, cvx_status] = NNLS_weighted_L1_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            objective = objective + (lambda1 * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), 1)) + (lambda2 * norm(x(indices(i,1)+1:indices(i,2)+1), inf)) ;
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end


% NNLS weighted group infinity norm + L1 penalty
function [x, cvx_status] = NNLS_L1_weighted_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b) + (lambda1 * norm(x,1));
        for i=1:size(groupWeights,1)
           objective = objective + (lambda2 * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), inf));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end


% NNLS weighted group infinity norm + weighted L1 penalty
function [x, cvx_status] = NNLS_weighted_L1_weighted_group_inf(A, b, m, groupWeights, indices, lambda1, lambda2)
    cvx_begin
        variable x(m);
        objective = norm(A*x-b);
        for i=1:size(groupWeights,1)
            objective = objective + (lambda1 * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), 1)) + (lambda2 * groupWeights(i) * norm(x(indices(i,1)+1:indices(i,2)+1), inf));
        end
        minimize(objective);
        subject to
            x >= 0;
    cvx_end
end

% NNLS sparse group lasso
function [x, cvx_status] = NNLS_sparse_group_lasso(A, b, m, n, indices, lambda1, lambda2, alpha, precursorOptions, fudge)
    TIC = sum(b);
    tot_option_abundance = 0.0;
    
    for i=1:size(indices,1)
        tot_option_abundance = tot_option_abundance + precursorOptions(i, 8);
    end

    cvx_begin
        cvx_solver_settings('maxit', 1000);
     
        variable x(m);
        %objective = 0;
        objective = norm(A*x-b);
        for i=1:size(indices,1)
            
            i1 = indices(i,1)+1;
            i2 = indices(i,2)+1;
            xi = x(i1:i2);
            p = sqrt(i2-i1);
            
            %objective = objective + (1/(2*n))*sum_square_pos(norm(A(:,i1:i2)*x(i1:i2) - b)) +  groupWeights(i)*(1-alpha)*lambda*p*(norm(x(i1:i2))) + groupWeights(i)*alpha*lambda*norm(x(i1:i2),1);
            objective = objective + lambda1 * norm(A(:,i1:i2)*x(i1:i2) - b) + (1-alpha)*lambda2*(norm(x(i1:i2)) + alpha*lambda2*norm(x(i1:i2),1));
            %objective = objective + groupWeights(i)*(1-alpha)*lambda*p*(norm(x(i1:i2))) + groupWeights(i)*alpha*lambda*norm(xi,1);
            
        end
        minimize(objective);
        subject to
            x >= min(b(b>0));
            for i=1:size(indices,1)
                i1 = indices(i,1)+1; 
                i2 = indices(i,2)+1;
                xi = x(i1:i2);
                Ai = A(:,i1:i2);

                norm(Ai*xi) <= (precursorOptions(i,8)/tot_option_abundance) * TIC * fudge;
            end
    cvx_end

    sprintf("TIC: %f", sum(A*x)/sum(b)) 
end

% NNLS sparse group lasso2                                                                                                                                  
function [x, cvx_status] = NNLS_sparse_group_lasso2(A, b, m, n, indices, lambda1, lambda2, alpha, precursorOptions, fudge)
    TIC = sum(b);

    tot_option_abundance = 0.0;

    for i=1:size(indices,1)
        tot_option_abundance = tot_option_abundance + precursorOptions(i, 8);
    end

    %lambda1 = lambda1 / 2^(size(indices,1) - 2);
    %lambda2 = lambda2 / 2^(size(indices,1) - 2);
    %lambda1 = lambda1 / (size(indices,1) - 1);
    %lambda2 = lambda2 / (size(indices,1) - 1); 


    cvx_begin
        variable x(m);
        %objective = (1-alpha)*lambda2*norm(x,2) + alpha*lambda2*norm(x,1);
        objective = norm(A*x-b) + lambda2*norm(x,1);     % + (1-alpha)*lambda2*norm(x,2) + alpha*lambda2*norm(x,1);
        %objective = lambda2*norm(x,1);



%        for i=1:m
%            posA = A(:,i) > 0;
%            for k=1:size(posA, 1)
%                j = posA(k);
%                if b(j) > 0
%                    objective = objective + lambda1 * norm(A(j,i) * x(i) - b(j), 2);
%                end
%            end
%        end


        for i=1:m
            posA = A(:,i) > 0;            
            objective = objective + lambda1 * norm(A(posA,i)*x(i) - b(posA));
        end




        %for i=1:size(indices,1)
        %     i1 = indices(i,1)+1;
        %     i2 = indices(i,2)+1;
        %     xi = x(i1:i2);
        %     
        %     objective = objective + 0.01 * norm(xi, inf);
        %end

        minimize(objective);
        subject to
            x >= 0; 
            %for i=1:size(indices,1) 
            %    i1 = indices(i,1)+1;
            %    i2 = indices(i,2)+1;
            %    xi = x(i1:i2);
            %    Ai = A(:,i1:i2);
            %    
            %    norm(Ai*xi) <= (precursorOptions(i,8)/tot_option_abundance) * TIC * fudge;
            % end
    cvx_end

    sprintf("TIC: %f", sum(A*x)/sum(b))     
end

% NNLS with relative abundance penalty
function [x, cvx_status] = NNLS_A(A, b, m, lambda1, precursorOptions, indices)
    cvx_begin
        variable x(m);

        minimize( norm(A*x-b) );


        subject to
            x >= 0;
    cvx_end
end

function writeRSpectrum(spectrum, mzValues, scanDetails, fileID, name)
    for j=1:size(spectrum,1)
        if spectrum(j) >= 1
            fprintf(fileID, '%f\t%f\t%d\t%s\n', mzValues(j)-0.0001, 0.0, scanDetails(1), name);
            fprintf(fileID, '%f\t%f\t%d\t%s\n', mzValues(j), spectrum(j), scanDetails(1), name);
            fprintf(fileID, '%f\t%f\t%d\t%s\n', mzValues(j)+0.0001, 0.0, scanDetails(1), name); 
        end
    end
end

function writeR(spectra, mzValues, precursorOptions, scanDetails, outFile, b)
    writeRSpectrum(b, mzValues, scanDetails, outFile, "observed");
    spectra(spectra<1) = 0;
    for i=1:size(spectra,2)
        if sum(spectra(:, i)) > min(b(b>0))
            charge = precursorOptions(i,5);
            
            title = strcat(['demixed=', num2str(i), ' charge=', num2str(charge), ' mass=', num2str(precursorOptions(i,1)), ' iso=', num2str(precursorOptions(i,3)), '-', num2str(precursorOptions(i,4))]);
            writeRSpectrum(spectra(:,i), mzValues, scanDetails, outFile, title);
        end
    end
end 


function writeMGF2(spectra, monoSpectra, mzValues, precursorOptions, calcPrecursorMass, scanDetails, outFile, globalTol, copy, b, minObs)
    spectra(spectra<minObs) = 0;
    monoSpectra(monoSpectra<minObs) = 0;
    for i=1:size(spectra,2)
        if sum(spectra(:, i)) > min(b(b>0))
            if copy
                writeMGFSpectrum(b, precursorOptions(i,:), mzValues, scanDetails, num2str(i), calcPrecursorMass, outFile, globalTol);
            else
                if precursorOptions(i,3) == 0
                    writeMGFSpectrum(spectra(:,i), precursorOptions(i,:), mzValues, scanDetails, num2str(i), calcPrecursorMass, outFile, globalTol);
                else
                    writeMGFSpectrum(monoSpectra(:,i), precursorOptions(i,:), mzValues, scanDetails, num2str(i), calcPrecursorMass, outFile, globalTol);
                end
            end
        end
    end
end

function writeMGF(spectra, mzValues, precursorOptions, calcPrecursorMass, scanDetails, outFile, globalTol, copy, b, minObs)
    spectra(spectra<1) = 0; 
    %spectra(spectra<minObs) = 0;
    for i=1:size(spectra,2)
        if sum(spectra(:, i)) > min(b(b>0))
            if copy
                writeMGFSpectrum(b, precursorOptions(i,:), mzValues, scanDetails, num2str(i), calcPrecursorMass, outFile, globalTol);  
            else
                writeMGFSpectrum(spectra(:,i), precursorOptions(i,:), mzValues, scanDetails, num2str(i), calcPrecursorMass, outFile, globalTol);    
            end
        end
    end
end

function writeOriginalMGF(b, mzValues, scanDetails, outFile, globalTol)
    intensity = 0;
    rt = scanDetails(4);
    charge = scanDetails(3);
    mass = scanDetails(2);
    title = strcat(['scan=', num2str(scanDetails(1)), ' demixed=original, charge=', num2str(charge), ' mass=', num2str(mass)]);
    nativeID = title;
    tol = num2str(globalTol); 

    filename = strcat(num2str(scanDetails(1)), '_ori_', tol, '.mgf');

    writeMGFHeader(outFile, title, filename, nativeID, mass, intensity, charge, rt, b, mzValues, scanDetails(1), 0);
end

function writeMGFHeader(fileID, title, filename, nativeID, mass, intensity, charge, rt, spectrum, mzValues, scanID, demixID)
    fprintf(fileID, 'BEGIN IONS\n');
    fprintf(fileID, 'TITLE=%s File:\"%s\", NativeID:\"%s\"\n', title, filename, nativeID);
    %fprintf(fileID, 'RTINSECONDS=%s\n', rt);
    fprintf(fileID, 'PEPMASS=%f %f\n', mass, intensity);
    fprintf(fileID, 'CHARGE=%d+\n', charge);
    fprintf(fileID, 'SCANS=%d\n', scanID*100 + demixID);

    for j=1:size(spectrum,1)
        if spectrum(j) >= 1
            fprintf(fileID, '%f %f\n', mzValues(j), spectrum(j));
        end
    end

    fprintf(fileID, 'END IONS\n');
end

function writeMGFSpectrum(spectrum, precursorOption, mzValues, scanDetails, i, calcPrecursorMass, outFile, globalTol)

    intensity = 0;
    rt = scanDetails(4);

    if calcPrecursorMass
        charge = precursorOption(5);
        mass = precursorOption(7);
        tol = num2str(globalTol+(precursorOption(2)-precursorOption(1))/2);
        title = strcat(['scan=', num2str(scanDetails(1)), ' demixed=', i, ' charge=', num2str(charge), ' minMass=', num2str(precursorOption(1)), ' maxMass=', num2str(precursorOption(2)), ' minIso=', num2str(precursorOption(3)), ' maxIso=', num2str(precursorOption(4))]);
        nativeID = title;
    else
        charge = scanDetails(3);
        mass = scanDetails(2);
        tol = num2str(globalTol);
        title = strcat(['scan=', num2str(scanDetails(1)), ' demixed=', i, ' charge=', num2str(precursorOption(5)), ' minMass=', num2str(precursorOption(1)), ' maxMass=', num2str(precursorOption(2)), ' minIso=', num2str(precursorOption(3)), ' maxIso=', num2str(precursorOption(4))]);
        nativeID = strcat(['scan=', num2str(scanDetails(1)), ' demixed=', i, ' charge=', num2str(charge), ' mass=', num2str(mass)]);
    end

    filename = strcat(num2str(scanDetails(1)), '_', i, '_', tol, '.mgf');

    writeMGFHeader(outFile, title, filename, nativeID, mass, intensity, charge, rt, spectrum, mzValues, scanDetails(1), str2num(i));

end


function spectra = computeSpectra(A, x, indices)
    spectra = zeros(size(A,1), size(indices,1));
    for i=1:size(indices,1)
        spectra(:,i) = A(:, indices(i,1)+1:indices(i,2)+1) * x(indices(i,1)+1:indices(i,2)+1);
    end    
end

function spectra = computeMonoSpectra(A, x, indices)
    spectra = zeros(size(A,1), size(indices,1));
    for i=1:size(indices,1)
        Ai =  A(:, indices(i,1)+1:indices(i,2)+1);
        xi = x(indices(i,1)+1:indices(i,2)+1);
        spectra(:,i) = computeBinnedMonoSpectrum(Ai, xi);
    end
end

function b = computeBinnedMonoSpectrum(A, x)
    [n, m] = size(A);
    b = zeros(n, 1);
    
    for i=1:m
        monoIndex = getMinNonNegativeIndex(A(:,i));
        b(monoIndex) = b(monoIndex) + norm(A(:,i), 1) * x(i);
    end
    
    b(b<1) = 0;
end

function i = getMinNonNegativeIndex(A)
    i = find(A>0, 1, 'first');
end


function outputCoefficients(A, x, b, indices, precursorOptions, algName, expName)
      
    sum(A*x)/sum(b)

    set(gcf, 'PaperPosition', [0 0 20 200])    % can be bigger than screen 
    set(gcf, 'PaperSize', [20 200])    % Same, but for PDF output
    
    f = figure('visible','off','Unit','inches','pos',[0 0 20 200]);
    % plot chimeric spectrum
    TIC = sum(b);
    subplot(size(indices,1)+1, 1, 1)
    bar(b)  
    set(gca, 'XTick', []) % remove X-axis ticks
    title({[expName, ' ', algName];
        ['TIC: ', num2str(TIC, '%g')]}, 'Interpreter', 'none', 'FontSize', 16);
    
    % evaluate each spectra
    for i=1:size(indices,1)
        spectrum = A(:, indices(i,1)+1:indices(i,2)+1) * x(indices(i,1)+1:indices(i,2)+1);
        % replace values <1 with 0
        spectrum(spectrum<0) = 0;
        TIC = sum(spectrum);
        maxTemplate = max(x(indices(i,1)+1:indices(i,2)+1));
        basePeak = max(spectrum);
        
        subplot(size(indices,1)+1, 1, i+1)
        bar(spectrum, 'black')   
        set(gca, 'XTick', []) % remove X-axis ticks
        title(opt2title(precursorOptions(i,:), TIC, basePeak, maxTemplate), 'FontSize', 16);
    end

    saveas(f, strcat('demixed_', expName, '_', algName), 'epsc')
end


function plotSpectra(spectra, b, precursorOptions, algName, expName, path)
    set(gcf, 'PaperPosition', [0 0 20 200])    % can be bigger than screen 
    set(gcf, 'PaperSize', [20 200])    % Same, but for PDF output
    
    f = figure('visible','off','Unit','inches','pos',[0 0 20 200]);
    % plot chimeric spectrum
    TIC = sum(b);
    subplot(size(spectra,2)+1, 1, 1)
    bar(b)  
    set(gca, 'XTick', []) % remove X-axis ticks
    title({[expName, ' ', algName];
        ['TIC: ', num2str(TIC, '%g')]}, 'Interpreter', 'none', 'FontSize', 16);
    
    % evaluate each spectra
    for i=1:size(spectra,2)
        spectrum = spectra(:,i);
        % replace values <1 with 0
        spectrum(spectrum<1) = 0;
        TIC = sum(spectrum);
        basePeak = max(spectrum);
        
        subplot(size(spectra,2)+1, 1, i+1)
        bar(spectrum, 'black')   
        set(gca, 'XTick', []) % remove X-axis ticks
        title(opt2title(precursorOptions(i,:), TIC, basePeak, 0), 'FontSize', 16);
    end

    saveas(f, strcat(path,'/demixed_', expName, '_', algName), 'epsc')
end


function x = opt2title(opt, TIC, basePeak, maxTemplate)
    x = {['Mass: ', num2str(opt(1)), ' - ', num2str(opt(2))]; 
         ['Charge: ', num2str(opt(5))];
         ['Isotopes: ', num2str(opt(3)), ' - ', num2str(opt(4))];
         ['Likelihood: ', num2str(opt(6)), ' Max Template: ', num2str(maxTemplate, '%g')]; 
         ['TIC: ', num2str(TIC, '%g'), ' BasePeak: ', num2str(basePeak, '%g')]};
end
