function []=Subplot_performance(params)

    fig = figure(1);
    fig.Position = [0, 0, 700, 600];

    % Define colors and markers for each method
    colors = {'b+-', 'R+-', 'k+-', 'c+-', 'yo-'};
    labels = {'Method2', 'My-Method2', 'My-Method2-ratio', 'Method2-111'};

    % Field names for M1_x, M2_x, ..., M5_x
    field_names = {'M1_', 'M2_', 'M3_', 'M4_', 'M5_'};
    y_labels = {'CD', 'GD', 'BD', 'ED'};  % Y-axis labels

    for i = 1:4
        subplot(2, 2, i);
        hold on;
        
        for j = 1:length(labels)
            field_name = [field_names{j}, num2str(i)];
            plot(params.SNR_dB, params.(field_name), colors{j});  
        end
        
        grid on;
        xlabel('SNR (dB)', 'Fontsize', 12);
        ylabel(y_labels{i}, 'Fontsize', 12);
        legend(labels, 'Location','southeast');
    end
    
end

