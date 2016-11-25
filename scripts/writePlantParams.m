% writes the plant parameters in a way CRootBox can read it

fname = 'test.pparam'; % choose the filename (file extension = pparam)

fid = fopen([pwd '/' fname],'w');

%m = monocot; % rename

fprintf(fid,'plantingdepth\t%g\n',m.plantingdepth);
fprintf(fid,'firstB\t%g\n',m.firstB);
fprintf(fid,'delayB\t%g\n',m.delayB);
fprintf(fid,'maxB\t%g\n',m.maxB);
fprintf(fid,'nC\t%g\n',m.nCR);
fprintf(fid,'firstSB\t%g\n',m.firstCR);
fprintf(fid,'delaySB\t%g\n',m.delayCR);
fprintf(fid,'delayRC\t%g\n',m.delayRC);
fprintf(fid,'nz\t%g\n',m.dzRC);

fclose(fid);
