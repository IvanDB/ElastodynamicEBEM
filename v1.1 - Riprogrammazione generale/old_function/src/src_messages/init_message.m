function fid = init_message(fid,domain_phys,domain_type,lev,module_type,...
    deg_k,kappa,velC,alpha,T_fin)

time_date=clock;

if strcmp(domain_phys,'time')
    
    fprintf('***   Starting TIGROU - Time-domain Problem   ***\n\n');
    fprintf('Domain: %s\n',domain_phys);
    fprintf('Mesh filename: %s\n',domain_type);
    fprintf('Level of refinement: %g\n',lev);
    fprintf('Speed of wave propagation: %g\n',velC);
    fprintf('Damping parameter: %g\n',alpha);
    fprintf('Final Time: %g\n',T_fin);
    fprintf('Method: %s\n',module_type);
    fprintf('Degree: %g\n',deg_k);
    
    fprintf('\n');
    
    fprintf('Date (Start): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
    fprintf('Time (Start): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
    
    fprintf(fid,'***   Starting TIGROU - Time-domain Problem   ***\n\n');
    fprintf(fid,'Domain: %s\n',domain_phys);
    fprintf(fid,'Mesh filename: %s\n',domain_type);
    fprintf(fid,'Level of refinement: %g\n',lev);
    fprintf(fid,'Speed of wave propagation: %g\n',velC);
    fprintf(fid,'Damping parameter: %g\n',alpha);
    fprintf(fid,'Final Time: %g\n',T_fin);
    fprintf(fid,'Method: %s\n',module_type);
    fprintf(fid,'Degree: %g\n',deg_k);
    
    fprintf(fid,'\n');
    
    fprintf(fid,'Date (Start): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
    fprintf(fid,'Time (Start): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
    
elseif strcmp(domain_phys,'freq')
    
    fprintf('***   Starting TIGROU - Frequency-domain Problem   ***\n\n');
    fprintf('Domain: %s\n',domain_phys);
    fprintf('Mesh filename: %s\n',domain_type);
    fprintf('Level of refinement: %g\n',lev);
    fprintf('Wave number: %g\n',kappa);
    fprintf('Method: %s\n',module_type);
    fprintf('Degree: %g\n',deg_k);
    
    fprintf('\n');
    
    fprintf('Date (Start): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
    fprintf('Time (Start): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
    
    fprintf(fid,'***   Starting TIGROU - Frequency-domain Problem   ***\n\n');
    fprintf(fid,'Domain: %s\n',domain_phys);
    fprintf(fid,'Mesh filename: %s\n',domain_type);
    fprintf(fid,'Level of refinement: %g\n',lev);
    fprintf(fid,'Wave number: %g\n',kappa);
    fprintf(fid,'Method: %s\n',module_type);
    fprintf(fid,'Degree: %g\n',deg_k);
    
    fprintf('\n');
    
    fprintf(fid,'Date (Start): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
    fprintf(fid,'Time (Start): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
    
end

end

