function [] = util_csvwrite(filename, x, str)

fid = fopen(filename,'wt'); % 書き込み用にファイルオープン
[rows,cols] = size(x);
if cols > 1
  fprintf(fid, '%s,', str{1,1:end-1}); % 文字列の書き出し
  fprintf(fid, '%s\n', str{1,end}); % 行末の文字列は、改行を含めて出力
  for i = 1:rows
      fprintf(fid, '%f,', x(i,1:end-1));
      fprintf(fid, '%f\n', x(i, end));
  end
else
  fprintf(fid, '%s\n', str{1});
  for i = 1:rows
      fprintf(fid, '%f\n', x(i));
  end
end
fclose('all');
 % ファイルクローズ
