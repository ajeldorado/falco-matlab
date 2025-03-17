function plot_n(materials,wvls)
   % plot_n(materials, wvls)
   % returns a plot of the refractive index of specific material(s)
   % accepts list of materials (list of str) or "all"
   
   if class(materials)=="char"
       disp("Use string as input (double quotes), not char (single quotes).")
       return
   end

   if materials=="all"
       mats = dir("*.csv");  % use all *.csv names
       materials=[];
       for i=1:length(mats)
           materials = [materials,string(mats(i).name(1:end-4))];
       end
   end
   
   figure();
   ax = gca;
   ax.FontSize = 15;
   ax.LineWidth = 1;
   hold on
   t = "Refractive index of: ";
   for material = materials
       n = eval(material+"(wvls)");
       plot(wvls,n,'LineWidth',2)
       t = t+material+", ";
   end
   xlabel('wavelength [µm]');
   ylabel('refractive index n');
   xlim([min(wvls),max(wvls)]);
   title(t)
   legend(materials)
end