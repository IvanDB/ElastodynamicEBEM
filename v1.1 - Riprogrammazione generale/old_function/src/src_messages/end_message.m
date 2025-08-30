function end_message(fid)

time_date=clock;


fprintf('Date (End): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
fprintf('Time (End): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
fprintf(fid,'\n');
fprintf('*** TIGROU has ended ***\n');

fprintf(fid,'Date (End): %d-%d-%d\n',time_date(3),time_date(2),time_date(1));
fprintf(fid, 'Time (End): %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));
fprintf(fid,'\n');
fprintf(fid,'*** TIGROU has ended ***\n');

end

