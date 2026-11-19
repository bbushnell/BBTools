package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Constructs symmetric integer matrices from retained directed monomer counts.
 * Pools counts before normalizing; never derives high-resolution scores from rounded scores.
 * @author Sayu
 */
public final class FamilySubstitutionMatrixExport{

	private FamilySubstitutionMatrixExport(){}
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args,FamilySubstitutionMatrixExport.class,false);
		try{run(pp.args);}finally{Shared.closeStream(pp.outstream);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg,'=');
			require(pair.length==2 && !pair[1].isEmpty(),"Require nonempty flag=value arguments");
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require((key.equals("in") || key.equals("insha80") || key.equals("out") || key.equals("scale") || key.equals("t")) &&
				options.put(key,pair[1])==null,"Unknown or duplicate matrix export option: "+key);
		}
		require(options.containsKey("in") && options.containsKey("insha80") && options.containsKey("out"),"Require in, insha80 and out");
		require(!options.containsKey("scale") || options.get("scale").equals("4"),"This recipe uses fourfold BLOSUM62 score precision");
		require(!options.containsKey("t") || options.get("t").equals("1"),"The 1200-row aggregate export uses one worker");
		final Path input=Paths.get(options.get("in")).toAbsolutePath().normalize();
		final Path output=Paths.get(options.get("out")).toAbsolutePath().normalize();
		final String pin=options.get("insha80");DigestSuffix.requireSuffix(pin,"Measured substitution counts");
		require(Files.isRegularFile(input) && !Files.isSymbolicLink(input),"Counts must be a regular input file");
		require(!Files.exists(output) && Files.isDirectory(output.getParent()),"Matrix output must be a fresh directory under an existing parent");
		require(pin.equals(DigestSuffix.file(input.toString())),"Measured count table differs from its declared sha80");
		final Counts counts=read(input);
		final Matrix low=construct(counts,1),pooled=construct(counts,2);
		require(pin.equals(DigestSuffix.file(input.toString())),"Measured count table changed during construction");
		Files.createDirectory(output);
		write(output.resolve("custom40.tsv"),matrixText("custom40",low));
		write(output.resolve("custom60.tsv"),matrixText("custom60",pooled));
		ProteinSubstitutionMatrix.load("custom40",output.resolve("custom40.tsv").toString());
		ProteinSubstitutionMatrix.load("custom60",output.resolve("custom60.tsv").toString());
		final ByteBuilder cells=new ByteBuilder("matrix\taa1\taa2\tunordered_count\tpairs\tendpoint1\tendpoint2\tobserved\texpected\tlog_odds_bits\tscore\n");
		appendCells(cells,"custom40",low);appendCells(cells,"custom60",pooled);
		write(output.resolve("cells.tsv"),cells);
		final ByteBuilder manifest=new ByteBuilder("key\tvalue\nsource\t").append(input.toString()).nl()
			.append("source_sha80\t").append(pin).nl().append("weighting\traw_eligible_columns\n")
			.append("custom40_bands\t<40\ncustom60_bands\t<40,40-60\n")
			.append("custom40_pairs\t").append(low.total).nl().append("custom60_pairs\t").append(pooled.total).nl()
			.append("custom40_alignments\t").append(low.alignments).nl().append("custom60_alignments\t").append(pooled.alignments).nl()
			.append("formula\tround_half_away(8*log2(observed/expected))\nscale\t4\ngap\t16\n")
			.append("smoothing\tnone;zero_cells_fail\nclipping\tnone\nnonstandard\tBLOSUM62_X_times4\n")
			.append("custom40_sha80\t").append(DigestSuffix.file(output.resolve("custom40.tsv").toString())).nl()
			.append("custom60_sha80\t").append(DigestSuffix.file(output.resolve("custom60.tsv").toString())).nl();
		write(output.resolve("manifest.tsv"),manifest);
		write(output.resolve("PASS"),new ByteBuilder("FAMILY_SUBSTITUTION_MATRIX_EXPORT_PASS matrices=2 unique_cells=420\n"));
	}

	/** Validate the full directed table and every repeated denominator before pooling. */
	private static Counts read(Path path){
		final Counts c=new Counts();final LineParser1 lp=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile(path.toString(),false);int rows=0;
		try{
			final byte[] header=reader.nextLine();
			require(header!=null && new String(header,StandardCharsets.US_ASCII).equals(HEADER),"Expected exact directed-count schema");
			for(byte[] line=reader.nextLine();line!=null;line=reader.nextLine()){
				lp.set(line);require(lp.terms()==8,"Directed-count row requires eight fields");
				final int band=lp.termEquals("<40",0) ? 0 : lp.termEquals("40-60",0) ? 1 : lp.termEquals(">=60",0) ? 2 : -1;
				require(band>=0,"Unexpected identity band");
				final int i=residue(lp.parseString(1)),j=residue(lp.parseString(2));
				require(!c.seen[band][i][j],"Duplicate directed cell in band "+band);c.seen[band][i][j]=true;
				final long n=Long.parseLong(lp.parseString(3)),q=Long.parseLong(lp.parseString(4)),r=Long.parseLong(lp.parseString(5)),
					a=Long.parseLong(lp.parseString(6)),t=Long.parseLong(lp.parseString(7));
				require(n>=0 && q>=0 && r>=0 && a>0 && t>0,"Counts must be nonnegative with positive population denominators");
				if(c.rows[band]==0){c.total[band]=t;c.alignments[band]=a;}
				require(c.total[band]==t && c.alignments[band]==a,"Repeated denominators disagree within a band");
				c.count[band][i][j]=n;c.query[band][i][j]=q;c.consensus[band][i][j]=r;c.rows[band]++;rows++;
			}
		}finally{require(!reader.close(),"Failed reading directed substitution counts");}
		require(rows==1200,"Each of three bands must contain all 400 directed cells");
		for(int b=0;b<3;b++){
			require(c.rows[b]==400,"Incomplete directed band: "+b);long total=0;
			for(int i=0;i<20;i++){
				long q=0,r=0;
				for(int j=0;j<20;j++){q=Math.addExact(q,c.count[b][i][j]);r=Math.addExact(r,c.count[b][j][i]);}
				for(int j=0;j<20;j++){
					require(c.query[b][i][j]==q && c.consensus[b][j][i]==r,"Repeated monomer counts must equal directed row/column sums");
				}
				total=Math.addExact(total,q);
			}
			require(total==c.total[b],"Directed counts do not conserve the eligible-pair denominator");
		}
		return c;
	}

	/** Pool directed counts first, then fold off-diagonals and derive endpoint-frequency log odds. */
	private static Matrix construct(Counts source,int bands){
		assert(bands==1 || bands==2) : "Only strict below40 and pooled below60 are released";
		final Matrix m=new Matrix();
		for(int b=0;b<bands;b++){
			m.total=Math.addExact(m.total,source.total[b]);m.alignments=Math.addExact(m.alignments,source.alignments[b]);
			for(int i=0;i<20;i++){for(int j=0;j<20;j++){m.count[i][j]=Math.addExact(m.count[i][j],source.count[b][i][j]);}}
		}
		for(int i=0;i<20;i++){
			for(int j=0;j<20;j++){m.endpoints[i]=Math.addExact(m.endpoints[i],Math.addExact(m.count[i][j],m.count[j][i]));}
			require(m.endpoints[i]>0,"Every standard residue needs positive endpoint mass");
		}
		long sum=0;for(long n:m.endpoints){sum=Math.addExact(sum,n);}
		require(sum==Math.multiplyExact(2,m.total),"Symmetric marginals count each observed pair twice");
		for(int i=0;i<20;i++){for(int j=i;j<20;j++){
			final long u=i==j ? m.count[i][j] : Math.addExact(m.count[i][j],m.count[j][i]);
			require(u>0,"Zero observed cell needs an explicit smoothing decision: "+ALPHABET.charAt(i)+ALPHABET.charAt(j));
			final double observed=u/(double)m.total,expected=(m.endpoints[i]/(2.0*m.total))*(m.endpoints[j]/(2.0*m.total))*(i==j ? 1 : 2);
			final double bits=StrictMath.log(observed/expected)/StrictMath.log(2);
			final int score=roundHalfAway(8*bits);
			require(i!=j || score>0,"Self-score normalization requires positive standard-residue diagonal scores");
			m.score[i][j]=m.score[j][i]=score;m.bits[i][j]=bits;m.observed[i][j]=observed;m.expected[i][j]=expected;
		}}
		return m;
	}

	/** Quantize signed values symmetrically, including exact positive and negative half ties. */
	public static int roundHalfAway(double value){
		require(Double.isFinite(value) && Math.abs(value)<Integer.MAX_VALUE-.5,"Log-odds score must fit an integer");
		return (int)Math.copySign(StrictMath.floor(Math.abs(value)+.5),value);
	}
	/** Emit explicit residue letters so table and BBTools encoded alphabet orders cannot be confused. */
	private static ByteBuilder matrixText(String name,Matrix m){
		assert(m.total>0) : "Matrix output requires a validated nonempty population";
		final ByteBuilder out=new ByteBuilder("#bbgap_substitution_matrix\t1\nname\t").append(name).append("\nscale\t4\ngap\t16\n");
		for(int i=0;i<20;i++){for(int j=0;j<20;j++){
			out.append("score\t").append(ALPHABET.charAt(i)).tab().append(ALPHABET.charAt(j)).tab().append(m.score[i][j]).nl();
		}}
		return out;
	}
	/** Retain each unordered count and its normalization terms beside the rounded score. */
	private static void appendCells(ByteBuilder out,String name,Matrix m){
		assert(m.total>0) : "Cell evidence uses the same validated denominator as matrix construction";
		for(int i=0;i<20;i++){for(int j=i;j<20;j++){
			out.append(name).tab().append(ALPHABET.charAt(i)).tab().append(ALPHABET.charAt(j)).tab()
				.append(i==j ? m.count[i][j] : Math.addExact(m.count[i][j],m.count[j][i])).tab().append(m.total).tab()
				.append(m.endpoints[i]).tab().append(m.endpoints[j]).tab().append(Double.toString(m.observed[i][j])).tab()
				.append(Double.toString(m.expected[i][j])).tab().append(Double.toString(m.bits[i][j])).tab().append(m.score[i][j]).nl();
		}}
	}
	private static int residue(String s){
		final int index=s.length()==1 ? ALPHABET.indexOf(s.charAt(0)) : -1;
		require(index>=0,"Only the 20 standard amino acids are permitted: "+s);return index;
	}
	private static void write(Path path,ByteBuilder data){
		require(!Files.exists(path),"Refusing to overwrite matrix evidence: "+path);
		final ByteStreamWriter out=new ByteStreamWriter(path.toString(),false,false,false);out.start();out.print(data);
		require(!out.poisonAndWait(),"Matrix output writer failed: "+path);
	}
	private static void require(boolean condition,String message){if(!condition){throw new IllegalArgumentException(message);}}
	private static final class Counts{
		final long[][][] count=new long[3][20][20],query=new long[3][20][20],consensus=new long[3][20][20];
		final boolean[][][] seen=new boolean[3][20][20];final long[] total=new long[3],alignments=new long[3];final int[] rows=new int[3];
	}
	private static final class Matrix{
		final long[][] count=new long[20][20];final long[] endpoints=new long[20];final int[][] score=new int[20][20];
		final double[][] bits=new double[20][20],observed=new double[20][20],expected=new double[20][20];long total,alignments;
	}
	private static final String ALPHABET="ARNDCEQGHILKMFPSTWYV";
	private static final String HEADER="band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs";
}
