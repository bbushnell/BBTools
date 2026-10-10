package prot;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.Parser;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Row-size analyzer for the experimental mincount HBM text dumps.
 *
 * @author Yelan
 */
public final class HbmTextDumpMincountAnalyzer {

	private HbmTextDumpMincountAnalyzer(){}

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest=t")){
			selfTest();
			return;
		}
		Path input=null, output=null;
		for(final String arg : Parser.parseConfig(args)){
			final int equals=arg.indexOf('=');
			if(equals<1 || equals==arg.length()-1){throw new IllegalArgumentException("Expected flag=value: "+arg);}
			final String key=arg.substring(0, equals), value=arg.substring(equals+1);
			if(key.equalsIgnoreCase("in") && input==null){input=Paths.get(value);}
			else if(key.equalsIgnoreCase("out") && output==null){output=Paths.get(value);}
			else{throw new IllegalArgumentException("Unknown or duplicate option: "+key);}
		}
		if(input==null || output==null){
			throw new IllegalArgumentException("Required in=<dump-directory> out=<fresh-output-directory>");
		}
		run(input, output);
	}

	static void run(final Path input, final Path output) throws Exception{
		if(!Files.isDirectory(input)){throw new IOException("Input directory not found: "+input);}
		if(Files.exists(output)){throw new IOException("Output already exists: "+output);}
		Files.createDirectories(output);
		final Analysis analysis=new Analysis();
		for(int minCount=1; minCount<=3; minCount++){
			analysis.analyzePaired(minCount, input.resolve("sparse_m"+minCount+".txt"),
				input.resolve("dense_m"+minCount+".txt"));
			analysis.analyzeImplicit(minCount, input.resolve("sparse_m"+minCount+".txt"));
		}
		analysis.write(output);
	}

	private static void selfTest() throws Exception{
		final Path root=Files.createTempDirectory("hbm-text-mincount-analyzer-");
		final Path input=root.resolve("input");
		final Path output=root.resolve("output");
		Files.createDirectories(input);
		final String denseZero=denseCounts(0, 0);
		for(int minCount=1; minCount<=3; minCount++){
			Files.write(input.resolve("sparse_m"+minCount+".txt"), sparseFixture().getBytes(StandardCharsets.UTF_8));
			Files.write(input.resolve("dense_m"+minCount+".txt"),
				denseFixture(denseZero).getBytes(StandardCharsets.UTF_8));
		}
		run(input, output);
		final String cardinality=read(output.resolve("hbm_text_mincount_row_cardinality_20261001.tsv"));
		final String summary=read(output.resolve("hbm_text_mincount_row_summary_20261001.tsv"));
		final String implicit=read(output.resolve("hbm_text_mincount_implicit_coordinate_estimate_20261001.tsv"));
		require(cardinality.contains("1\tREF\t2\t1\t14\t52\t6\t44\n"), "REF cardinality row wrong");
		require(cardinality.contains("1\tREF_INS\t0\t1\t10\t54\t0\t44\n"), "zero-state REF_INS row wrong");
		require(summary.contains("1\tALL\t2\t1\t24\t106\t6\t88\t18\t18\n"), "summary row wrong");
		require(implicit.contains("1\tDEL\t1\t6\t4\t2\t6\t4\n"), "implicit DEL row wrong");
		System.out.println("HbmTextDumpMincountAnalyzer PASS");
		System.out.println("fixture_directory\t"+root);
	}

	private static String sparseFixture(){
		return "#format\thbm_text_v1\n"+
			"#encoding\tsparse\n"+
			"#min_count_emitted\t1\n"+
			"f\t0\tfam\t1\tabc\tA\n"+
			"r\t0\tr\t5\tA1\tR2\n"+
			"i\t0\t0\tr\t2\n"+
			"d\t0\t:\n"+
			"e\n"+
			"z\n";
	}

	private static String denseFixture(final String zeroCounts){
		return "#format\thbm_text_v1\n"+
			"#encoding\tdense\n"+
			"#min_count_emitted\t1\n"+
			"f\t0\tfam\t1\tabc\tA\n"+
			"r\t0\td\t5"+denseCounts(1, 2)+"\n"+
			"i\t0\t0\td\t2"+zeroCounts+"\n"+
			"d\t0\t:\n"+
			"e\n"+
			"z\n";
	}

	private static String denseCounts(final int first, final int second){
		final StringBuilder sb=new StringBuilder();
		sb.append('\t').append(first).append('\t').append(second);
		for(int i=2; i<HbmBundleFormat.NAA; i++){sb.append("\t0");}
		return sb.toString();
	}

	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
	}

	private static void require(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final class Analysis {
		final RowStats[][][] cardinality=new RowStats[3][3][HbmBundleFormat.NAA+1];
		final RowStats[][] summary=new RowStats[3][3];
		final ImplicitStats[][] implicit=new ImplicitStats[3][5];
		final RowStats[] allSummary=new RowStats[3];
		final ImplicitStats[] allImplicit=new ImplicitStats[3];

		Analysis(){
			for(int min=0; min<3; min++){
				for(int i=0; i<cardinality[min].length; i++){
					for(int j=0; j<cardinality[min][i].length; j++){
						cardinality[min][i][j]=new RowStats();
					}
				}
				for(int i=0; i<summary[min].length; i++){summary[min][i]=new RowStats();}
				for(int i=0; i<implicit[min].length; i++){implicit[min][i]=new ImplicitStats();}
				allSummary[min]=new RowStats();
				allImplicit[min]=new ImplicitStats();
			}
		}

		void analyzePaired(final int minCount, final Path sparse, final Path dense) throws IOException{
			final ByteFile sbf=ByteFile.makeByteFile(sparse.toString(), false);
			final ByteFile dbf=ByteFile.makeByteFile(dense.toString(), false);
			final LineParser1 sp=new LineParser1('\t'), dp=new LineParser1('\t');
			try{
				byte[] sLine, dLine;
				long line=0;
				while((sLine=sbf.nextLine())!=null){
					line++;
					dLine=dbf.nextLine();
					if(dLine==null){throw new IOException("Dense file ended before sparse at line "+line);}
					analyzePairedLine(minCount, sLine, dLine, sp, dp, line);
				}
				if(dbf.nextLine()!=null){throw new IOException("Sparse file ended before dense");}
			}finally{
				if(sbf.close()){throw new IOException("Error closing "+sparse);}
				if(dbf.close()){throw new IOException("Error closing "+dense);}
			}
		}

		private void analyzePairedLine(final int minCount, final byte[] sLine, final byte[] dLine,
				final LineParser1 s, final LineParser1 d, final long line) throws IOException{
			s.set(sLine);
			d.set(dLine);
			if(sLine[0]=='#'){return;}
			if(s.termEquals('f', 0) || s.termEquals('d', 0) ||
					s.termEquals('e', 0) || s.termEquals('z', 0)){
				if(!Arrays.equals(sLine, dLine)){throw new IOException("Metadata mismatch at line "+line);}
				return;
			}
			final int rowType, fixed;
			if(s.termEquals('r', 0)){rowType=0; fixed=4;}
			else if(s.termEquals('i', 0)){rowType=1; fixed=5;}
			else if(s.termEquals("di", 0)){rowType=2; fixed=5;}
			else{throw new IOException("Unknown row tag at line "+line);}
			if(!sameTerm(s, d, 0)){throw new IOException("Row tag mismatch at line "+line);}
			assertShared(s, d, fixed, line);
			final int nonzero=s.terms()-fixed;
			if(nonzero<0 || nonzero>HbmBundleFormat.NAA){
				throw new IOException("Sparse cardinality out of range at line "+line+": "+nonzero);
			}
			final long sparseBytes=rowBytes(sLine), denseBytes=rowBytes(dLine);
			final long sparseState=stateBytes(s, fixed), denseState=stateBytes(d, fixed);
			final int minIndex=minCount-1;
			cardinality[minIndex][rowType][nonzero].add(nonzero, sparseBytes,
				denseBytes, sparseState, denseState);
			summary[minIndex][rowType].add(nonzero, sparseBytes, denseBytes, sparseState, denseState);
			allSummary[minIndex].add(nonzero, sparseBytes, denseBytes, sparseState, denseState);
		}

		void analyzeImplicit(final int minCount, final Path sparse) throws IOException{
			final ByteFile bf=ByteFile.makeByteFile(sparse.toString(), false);
			final LineParser1 lp=new LineParser1('\t');
			try{
				byte[] line;
				long lineNumber=0;
				while((line=bf.nextLine())!=null){
					lineNumber++;
					analyzeImplicitLine(minCount, line, lp, lineNumber);
				}
			}finally{
				if(bf.close()){throw new IOException("Error closing "+sparse);}
			}
		}

		private void analyzeImplicitLine(final int minCount, final byte[] line, final LineParser1 fields,
				final long lineNumber) throws IOException{
			fields.set(line);
			final int rowType;
			final long currentBytes=rowBytes(line), implicitBytes, stateBytes;
			if(line[0]=='#' || fields.termEquals('f', 0) ||
					fields.termEquals('e', 0) || fields.termEquals('z', 0)){
				rowType=0;
				implicitBytes=currentBytes;
				stateBytes=0;
			}else if(fields.termEquals('d', 0)){
				rowType=3;
				implicitBytes=3+fields.length(2);
				stateBytes=0;
			}else if(fields.termEquals('r', 0)){
				rowType=1;
				stateBytes=stateBytes(fields, 4);
				implicitBytes=stateBytes+5+fields.length(3);
			}else if(fields.termEquals('i', 0)){
				rowType=2;
				stateBytes=stateBytes(fields, 5);
				implicitBytes=stateBytes+5+fields.length(4);
			}else if(fields.termEquals("di", 0)){
				rowType=4;
				stateBytes=stateBytes(fields, 5);
				implicitBytes=stateBytes+6+fields.length(4);
			}else{
				throw new IOException("Unknown row tag at line "+lineNumber);
			}
			final int minIndex=minCount-1;
			implicit[minIndex][rowType].add(currentBytes, implicitBytes, stateBytes);
			allImplicit[minIndex].add(currentBytes, implicitBytes, stateBytes);
		}

		private static void assertShared(final LineParser1 s, final LineParser1 d,
				final int fixed, final long line) throws IOException{
			if(d.terms()!=fixed+HbmBundleFormat.NAA){
				throw new IOException("Dense count field mismatch at line "+line);
			}
			for(int i=0; i<fixed; i++){
				if(i==fixed-2){continue;}
				if(!sameTerm(s, d, i)){
					throw new IOException("Coordinate mismatch at line "+line);
				}
			}
			if(!s.termEquals('r', fixed-2) || !d.termEquals('d', fixed-2)){
				throw new IOException("Row encoding mismatch at line "+line);
			}
		}

		void write(final Path output) throws IOException{
			writeCardinality(output.resolve("hbm_text_mincount_row_cardinality_20261001.tsv"));
			writeSummary(output.resolve("hbm_text_mincount_row_summary_20261001.tsv"));
			writeImplicit(output.resolve("hbm_text_mincount_implicit_coordinate_estimate_20261001.tsv"));
		}

		private void writeCardinality(final Path path) throws IOException{
			final ByteStreamWriter bsw=writer(path);
			final ByteBuilder bb=new ByteBuilder(1<<16);
			try{
				bb.append("min_count\trow_type\tnonzero\trows\tsparse_row_bytes\tdense_row_bytes\tsparse_state_bytes\tdense_state_bytes").nl();
				for(int min=1; min<=3; min++){
					for(int type=0; type<ROW_TYPES.length; type++){
						for(int nonzero=0; nonzero<=HbmBundleFormat.NAA; nonzero++){
							final RowStats stats=cardinality[min-1][type][nonzero];
							if(stats.rows>0){
								bb.append(min).tab().append(ROW_TYPES[type]).tab();
								bb.append(nonzero).tab();
								stats.appendRow(bb).nl();
							}
						}
					}
				}
				bsw.print(bb);
			}finally{
				if(bsw.poisonAndWait()){throw new IOException("Error writing "+path);}
			}
		}

		private void writeSummary(final Path path) throws IOException{
			final ByteStreamWriter bsw=writer(path);
			final ByteBuilder bb=new ByteBuilder(1<<16);
			try{
				bb.append("min_count\trow_type\trows\tmean_nonzero\tsparse_row_bytes\tdense_row_bytes\t");
				bb.append("sparse_state_bytes\tdense_state_bytes\tsparse_coord_bytes\tdense_coord_bytes").nl();
				for(int min=1; min<=3; min++){
					for(int type=0; type<ROW_TYPES.length; type++){
						bb.append(min).tab().append(ROW_TYPES[type]).tab();
						summary[min-1][type].appendSummary(bb).nl();
					}
					bb.append(min).tab().append("ALL").tab();
					allSummary[min-1].appendSummary(bb).nl();
				}
				bsw.print(bb);
			}finally{
				if(bsw.poisonAndWait()){throw new IOException("Error writing "+path);}
			}
		}

		private void writeImplicit(final Path path) throws IOException{
			final ByteStreamWriter bsw=writer(path);
			final ByteBuilder bb=new ByteBuilder(1<<16);
			try{
				bb.append("min_count\trow_type\trows\tcurrent_bytes\timplicit_coordinate_bytes\t");
				bb.append("saved_bytes\tcurrent_coord_bytes\timplicit_coord_bytes").nl();
				for(int min=1; min<=3; min++){
					for(int type=0; type<IMPLICIT_TYPES.length; type++){
						bb.append(min).tab().append(IMPLICIT_TYPES[type]).tab();
						implicit[min-1][type].appendRow(bb).nl();
					}
					bb.append(min).tab().append("ALL").tab();
					allImplicit[min-1].appendRow(bb).nl();
				}
				bsw.print(bb);
			}finally{
				if(bsw.poisonAndWait()){throw new IOException("Error writing "+path);}
			}
		}
	}

	private static ByteStreamWriter writer(final Path path){
		final ByteStreamWriter bsw=new ByteStreamWriter(path.toString(), true, false, false);
		bsw.start();
		return bsw;
	}

	private static long rowBytes(final byte[] line){return line.length+1L;}

	private static long stateBytes(final LineParser1 fields, final int from){
		long bytes=0;
		for(int i=from; i<fields.terms(); i++){bytes+=1L+fields.length(i);}
		return bytes;
	}

	private static boolean sameTerm(final LineParser1 a, final LineParser1 b, final int term){
		final int length=a.length(term);
		if(b.length(term)!=length){return false;}
		for(int i=0; i<length; i++){
			if(a.parseByte(term, i)!=b.parseByte(term, i)){return false;}
		}
		return true;
	}

	private static String mean(final long sum, final long rows){
		if(rows==0){return "0";}
		final BigDecimal bd=BigDecimal.valueOf((double)sum/(double)rows)
			.round(new MathContext(6, RoundingMode.HALF_UP)).stripTrailingZeros();
		return bd.toPlainString();
	}

	private static final String[] ROW_TYPES={"REF", "REF_INS", "DEL_INS"};
	private static final String[] IMPLICIT_TYPES={"META", "REF", "REF_INS", "DEL", "DEL_INS"};

	private static final class RowStats {
		long rows;
		long nonzeroSum;
		long sparseRowBytes;
		long denseRowBytes;
		long sparseStateBytes;
		long denseStateBytes;

		void add(final int nonzero, final long sparseRow, final long denseRow,
				final long sparseState, final long denseState){
			rows++;
			nonzeroSum+=nonzero;
			sparseRowBytes+=sparseRow;
			denseRowBytes+=denseRow;
			sparseStateBytes+=sparseState;
			denseStateBytes+=denseState;
		}

		ByteBuilder appendRow(final ByteBuilder bb){
			return bb.append(rows).tab().append(sparseRowBytes).tab()
				.append(denseRowBytes).tab().append(sparseStateBytes).tab()
				.append(denseStateBytes);
		}

		ByteBuilder appendSummary(final ByteBuilder bb){
			return bb.append(rows).tab().append(mean(nonzeroSum, rows)).tab()
				.append(sparseRowBytes).tab().append(denseRowBytes).tab()
				.append(sparseStateBytes).tab().append(denseStateBytes).tab()
				.append(sparseRowBytes-sparseStateBytes).tab()
				.append(denseRowBytes-denseStateBytes);
		}
	}

	private static final class ImplicitStats {
		long rows;
		long currentBytes;
		long implicitBytes;
		long currentCoordBytes;
		long implicitCoordBytes;

		void add(final long current, final long implicit, final long state){
			rows++;
			currentBytes+=current;
			implicitBytes+=implicit;
			currentCoordBytes+=current-state;
			implicitCoordBytes+=implicit-state;
		}

		ByteBuilder appendRow(final ByteBuilder bb){
			return bb.append(rows).tab().append(currentBytes).tab()
				.append(implicitBytes).tab().append(currentBytes-implicitBytes).tab()
				.append(currentCoordBytes).tab().append(implicitCoordBytes);
		}
	}
}
