start_idx = 847;
end_idx = 848;

for idx = start_idx:end_idx
    fname = sprintf('%04d.npy', idx);   % Pads with zeros → '0847.npy'
    fprintf('\n===== Processing file: %s =====\n', fname);

    if ~isfile(fname)
        warning('File %s does not exist. Skipping...', fname);
        continue
    end

    img = readNPY(fname); 
    %imshow(squeeze(img(1, :, :)), []);
    
        for i = 1:size(img,1)
        I = squeeze(img(i,:,:));
        I = im2double(I); 
    
        % Start of the original code
        Nmodel = 256;
        Nacq = 256;
        readoutDirection = 'z';
        echoTime = 3;
        samplingInterval = 2;
        
        %Put your image's size the same as the model
        I = imresize(I, [Nmodel, Nmodel]);
        
        fprintf('========\n');
        fprintf('Example 2: Titanium cylinder in 3D\nModel matrix: %dx%dx%d\nAcquisition matrix: %dx%dx%d\nReadout direction: %s\nTE: %f ms\nSampling interval: %f ms\n', Nmodel, Nmodel, Nmodel, Nacq, Nacq, Nacq, readoutDirection, echoTime, samplingInterval);
        
        % carrotImage = createCarrotImage(Nmodel);
        
        % Create a cylindrical mask
        [x,y,z] = meshgrid(linspace(-1,1,Nmodel),linspace(-1,1,Nmodel),linspace(-1,1,Nmodel));
        
        % Define offsets (choose your own values)
        offsetX = 0.0;  % +0.2 moves cylinder in +x direction
        offsetY = 0.0; % -0.1 moves cylinder in -y direction
        offsetZ = -0.0;  % 0 => no shift in z
        
        % Now construct the mask with these offsets
        cylinderMask = sqrt((x - offsetX).^2 + (y - offsetY).^2) < 0.1 & abs(z - offsetZ) <= 0.3;
        
        % Create simulation model: Density 1 in background, density 0 inside cylinder
        model = struct();
        model.protonDensity = ones(Nmodel, Nmodel, Nmodel);
        for j = 1:Nmodel
            model.protonDensity(:,j,:) = I;
        end
        model.protonDensity = permute(model.protonDensity, [3 2 1]);
        model.protonDensity(cylinderMask) = 0;
        
        susWater = -9e-6;
        susTitanium = 180e-6;
        
        susceptibility = susWater * ones(Nmodel,Nmodel,Nmodel);
        susceptibility(cylinderMask) = susTitanium;
        
        % Calculate susceptibility relative to background susceptibility to avoid boundary artifacts
        susceptibility = susceptibility - susWater;
        
        model.resolution = [Nacq/Nmodel Nacq/Nmodel Nacq/Nmodel];
        fprintf('Calculating susceptibility-induced field shift... ');
        tic
        model.deltaB0 = calculateFieldShift(susceptibility, model.resolution);
        toc
        
        % Acquisition parameters
        acquisition = struct();
        acquisition.kspaceSamplingTimes = calculateCartesianSamplingTimes('gradientecho', [Nacq Nacq Nacq], readoutDirection, echoTime, samplingInterval);
        acquisition.resolution = [1 1 1];
        
        % Simulation
        fprintf('Simulating... ');
        tic
        kspace = forecast(model, acquisition);
        toc
        image = ifftc(kspace);
        image = 2 * image;
        
        
        figure
        subplot(2,2,1), imagesc(rot90(squeeze(model.protonDensity(:,fftCenter(Nmodel),:)),-1))
        colormap(gray(256))
        axis equal
        title('3D Titanium cylinder - Proton density')
        subplot(2,2,2), imagesc(rot90(squeeze(model.deltaB0(:,fftCenter(Nmodel),:)),-1))
        colormap(gray(256))
        axis equal
        title('Delta-B0')
        
        subplot(2,2,3), imagesc(rot90(abs(squeeze(image(:,fftCenter(Nacq),:))),-1),[0 1.5])
        colormap(gray(256))
        axis equal
        title('Z readout - Coronal slice - Magnitude')
        subplot(2,2,4), imagesc(rot90(angle(squeeze(image(:,fftCenter(Nacq),:))),-1),[-pi pi])
        colormap(gray(256))
        axis equal
        title('Phase')
        drawnow
    
        % Choose a central coronal slice for both outputs
        center_pd  = fftCenter(Nmodel);  % for proton density
        center_mag = fftCenter(Nacq);    % for magnitude FFT image
        
        % Extract 2D slices
        pd_slice   = squeeze(model.protonDensity(:, center_pd, :));      % size: [512 x 512]
        mag_slice  = squeeze(abs(image(:, center_mag, :)));              % size: [256 x 256]
        
        % Output path
        outdir = fullfile('outputs', sprintf('%04d', idx));
        if ~exist(outdir, 'dir'); mkdir(outdir); end
    
        outfile1 = fullfile(outdir, sprintf('protonDensity_slice%d.npy', i));
        outfile2 = fullfile(outdir, sprintf('magnitudeFFT_slice%d.npy', i));
    
        writeNPY(single(pd_slice), outfile1);
        writeNPY(single(mag_slice), outfile2);
        end
end