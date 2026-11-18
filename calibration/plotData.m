function plotData(data)
    figure('Position', [100 100 1200 400]);
    subplot(1, 2, 1);
    hold on;
    plot(data.time.all, data.annual.f, 'r-o', 'LineWidth', 2, 'DisplayName', 'Female');
    plot(data.time.all, data.annual.m, 'b-s', 'LineWidth', 2, 'DisplayName', 'Male');
    xlabel('Time (year)');
    ylabel('Number of Cases');
    title('Annual HPV-related Cancer Incidence');
    legend('Location', 'best');
    grid on;
    subplot(1, 2, 2);
    hold on;
    plot(data.time.all, data.annual.f_death, 'r-o', 'LineWidth', 2);
    xlabel('Time (year)');
    ylabel('Number of Deaths');
    title('Annual Cervical Cancer Mortality');
    ylim([0 1200]);
    grid on;
end
