package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Native and NCBI text interchange for canonical protein substitution scores.
 * NCBI ambiguity/stop cells are validated but never installed as runtime scores.
 * Export backgrounds are explicit probabilities in BBTools encoded residue order.
 * This class owns each opened reader/writer; all parsing state is invocation-local.
 * @author Fischl
 */
public final class ProteinSubstitutionMatrixIO{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static format and conversion operations only. */
	private ProteinSubstitutionMatrixIO(){}

	/** Converts one small matrix using ordinary BBTools argument and stream handling. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, ProteinSubstitutionMatrixIO.class, false);
		try{convert(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/** Detects content independently of extension. Null name and zero units mean
	 * unspecified; plain NCBI files require all three explicit values. Supplied
	 * values must agree with any native or NCBI metadata, never override it.
	 */
	public static ProteinSubstitutionMatrix read(String path, String name, int scale, int gap){
		require(path!=null && !path.isEmpty(), "Matrix resource path is required");
		require(scale>=0 && gap>=0, "Unspecified units use zero; negative units are invalid");
		final ByteFile input=ByteFile.makeByteFile1(path, false);
		try{
			final byte[] first=input.nextLine();
			if(starts(first, "#bbgap_substitution_matrix")){
				final ProteinSubstitutionMatrix matrix=ProteinSubstitutionMatrix.readNative(path, input, first);
				require(name==null || name.equals(matrix.name), "Conflicting matrix name metadata");
				require(scale==0 || scale==matrix.scale, "Conflicting matrix scale metadata");
				require(gap==0 || gap==matrix.linearGap, "Conflicting matrix gap metadata");
				return matrix;
			}
			return readNcbi(input, first, name, scale, gap);
		}finally{
			require(!input.close(), "Substitution matrix resource read failed: "+path);
		}
	}

	/** Parses labeled rows, permitting reordered axes and recognized extra symbols.
	 * Every declared row and column must occur once, even for ignored extras.
	 * The runtime factory separately checks canonical symmetry and arithmetic bounds.
	 */
	private static ProteinSubstitutionMatrix readNcbi(ByteFile input, byte[] first,
			String name, int scale, int gap){
		assert(input!=null) : "The format dispatcher owns this open byte reader";
		final LineParser1 parser=new LineParser1('\t');
		final ByteBuilder scratch=new ByteBuilder();
		final int[][] scores=new int[20][20];
		final boolean[] columns=new boolean[INPUT_SYMBOLS.length()], rows=new boolean[INPUT_SYMBOLS.length()];
		int[] order=null;
		int rowCount=0;
		boolean metadata=false;
		for(byte[] line=first; line!=null; line=input.nextLine()){
			final byte[] normalized=normalize(line, scratch);
			if(normalized.length==0){continue;}
			parser.set(normalized);
			if(normalized[0]=='#'){
				require(!parser.termStartsWith("#bbgap_substitution_matrix", 0),
					"Native matrix magic must be the exact first line with its supported version");
				if(parser.termStartsWith("#bbgap_matrix", 0)){
					require(!metadata && order==null && parser.terms()==5
						&& parser.termEquals("#bbgap_matrix", 0) && parser.termEquals("1", 1),
						"Invalid, repeated or misplaced NCBI matrix metadata version");
					require(parser.termStartsWith("name=", 2) && parser.termStartsWith("scale=", 3)
						&& parser.termStartsWith("gap=", 4), "NCBI metadata requires name, scale and gap");
					final String suppliedName=parser.parseString(2).substring(5);
					final int suppliedScale=integer(parser, 3, 6), suppliedGap=integer(parser, 4, 4);
					require(name==null || name.equals(suppliedName), "Conflicting matrix name metadata");
					require(scale==0 || scale==suppliedScale, "Conflicting matrix scale metadata");
					require(gap==0 || gap==suppliedGap, "Conflicting matrix gap metadata");
					name=suppliedName;
					scale=suppliedScale;
					gap=suppliedGap;
					metadata=true;
				}
				continue;
			}
			if(order==null){
				require(parser.terms()>=20 && parser.terms()<=INPUT_SYMBOLS.length(),
					"NCBI matrix header needs canonical20 and optional B/Z/X/star/J columns");
				order=new int[parser.terms()];
				for(int c=0; c<order.length; c++){
					final int symbol=symbol(parser, c);
					require(!columns[symbol], "Duplicate NCBI column label");
					columns[symbol]=true;
					order[c]=symbol;
				}
				for(int a=0; a<20; a++){require(columns[a], "Missing canonical NCBI column");}
			}else{
				require(parser.terms()==order.length+1, "NCBI row width differs from its header");
				final int a=symbol(parser, 0);
				require(columns[a] && !rows[a], "Undeclared or duplicate NCBI row label");
				rows[a]=true;
				rowCount++;
				for(int c=0; c<order.length; c++){
					final int value=integer(parser, c+1, 0), b=order[c];
					if(a<20 && b<20){scores[code(a)][code(b)]=value;}
				}
			}
		}
		require(order!=null && rowCount==order.length && Arrays.equals(rows, columns),
			"Incomplete NCBI matrix rows");
		require(name!=null && scale>0 && gap>0,
			"Plain NCBI input requires explicit name, scale and gap; units are never inferred from a filename");
		return ProteinSubstitutionMatrix.fromScores(name, scale, gap, scores);
	}

	/** Converts only ASCII whitespace separators to single tabs; leaves token bytes intact. */
	private static byte[] normalize(byte[] line, ByteBuilder scratch){
		assert(line!=null && scratch!=null) : "NCBI token normalization borrows a line and reusable scratch";
		scratch.clear();
		boolean separator=false;
		for(byte b:line){
			if(b==' ' || b=='\t' || b=='\r' || b=='\f' || b==11){separator=scratch.length()>0;}
			else{
				if(separator){scratch.tab();}
				scratch.append(b);
				separator=false;
			}
		}
		return scratch.toBytes();
	}

	/** Reads checked signed decimal integers without narrowing or overflow before validation. */
	private static int integer(LineParser1 parser, int term, int offset){
		parser.setBounds(term);
		final byte[] line=parser.line();
		int at=parser.a()+offset;
		final int end=parser.b();
		require(at<end, "Empty matrix integer");
		final boolean negative=line[at]=='-';
		if(negative || line[at]=='+'){at++;}
		require(at<end, "Matrix integer requires decimal digits");
		final long limit=negative ? 2147483648L : Integer.MAX_VALUE;
		long value=0;
		for(; at<end; at++){
			final int digit=line[at]-'0';
			require(digit>=0 && digit<=9, "Noninteger NCBI matrix cell or metadata");
			require(value<=(limit-digit)/10, "Matrix integer overflow");
			value=value*10+digit;
		}
		return (int)(negative ? -value : value);
	}

	/** Recognizes the canonical alphabet and five explicit interchange-only symbols. */
	private static int symbol(LineParser1 parser, int term){
		require(parser.length(term)==1, "NCBI labels must be single uppercase residue symbols");
		final int symbol=INPUT_SYMBOLS.indexOf((char)parser.line()[parser.a()]);
		require(symbol>=0, "Unknown NCBI residue label");
		return symbol;
	}

	/** Maps the external alphabet order to the established encoded canonical order. */
	private static int code(int index){
		assert(index>=0 && index<20) : "Only canonical columns map to runtime score coordinates";
		return AminoAcid.acidToNumberExtended[CANONICAL.charAt(index)];
	}

	/** Checks a small ASCII magic prefix without converting input to a String. */
	private static boolean starts(byte[] line, String prefix){
		if(line==null || line.length<prefix.length()){return false;}
		for(int i=0; i<prefix.length(); i++){if(line[i]!=prefix.charAt(i)){return false;}}
		return true;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Writing            ----------------*/
	/*--------------------------------------------------------------*/

	/** Emits native version1 with all400 canonical cells; no ambiguity cells or new required metadata. */
	public static ByteBuilder nativeText(ProteinSubstitutionMatrix matrix){
		require(matrix!=null, "Matrix required for native export");
		final ByteBuilder out=new ByteBuilder("#bbgap_substitution_matrix\t1\nname\t");
		out.append(matrix.name).append("\nscale\t").append(matrix.scale)
			.append("\ngap\t").append(matrix.linearGap).nl();
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){
				out.append("score").tab().append(CANONICAL.charAt(a)).tab().append(CANONICAL.charAt(b))
					.tab().append(matrix.score((byte)code(a), (byte)code(b))).nl();
			}
		}
		return out;
	}

	/** Emits canonical20+B/Z/X/star using an explicit background in encoded order.
	 * Ambiguous pairs are expectations over both residue distributions, rounded once
	 * half away from zero. Canonical scores are copied exactly. Background/source
	 * describe this export only and never alter the immutable runtime matrix.
	 */
	public static ByteBuilder ncbiText(ProteinSubstitutionMatrix matrix, double[] background, String source){
		return ncbiText(matrix, background, source, false);
	}

	/** Adds derivation claims only after the converter has matched the supplied
	 * probability state to every canonical score; a name alone never proves them.
	 */
	private static ByteBuilder ncbiText(ProteinSubstitutionMatrix matrix, double[] background,
			String source, boolean verifiedPossum){
		require(matrix!=null && source!=null && !source.isEmpty(), "Matrix and background source are required");
		for(int i=0; i<source.length(); i++){
			require(source.charAt(i)>=32 && source.charAt(i)<127, "Background source must be one printable ASCII line");
		}
		final double[][] weights=weights(background);
		final ByteBuilder out=new ByteBuilder("#bbgap_matrix\t1\tname=");
		out.append(matrix.name).append("\tscale=").append(matrix.scale).append("\tgap=").append(matrix.linearGap).nl();
		if(verifiedPossum){
			require(ProteinSubstitutionMatrix.POSSUM55_NAME.equals(matrix.name), "POSSUM55 description requires its name");
			require(matrix.scale==4 && matrix.linearGap==16, "POSSUM55 export requires native scale4 and gap16");
			out.append("# POSSUM55: ").append(ProteinSubstitutionMatrix.POSSUM55_DESCRIPTION).nl();
			out.append("# Identity band: 25-55%; scores round(8*log2(observed/expected)), eighth-bit units.\n");
			out.append("# Supplied25to55 probability state reproduces every canonical score; source provenance is declared below.\n");
		}
		out.append("# Linear gap penalty: ").append(matrix.linearGap).append(" score units.\n");
		out.append("# Background source: ").append(source).nl();
		out.append("# Background probabilities:");
		for(int a=0; a<20; a++){
			out.append(' ').append(CANONICAL.charAt(a)).append('=').appendSlow(background[code(a)]);
		}
		out.nl().append("# B=D/N, Z=E/Q: background-weighted; X: background expected score; round once, half away from zero.\n");
		out.append("# Star/residue: canonical matrix minimum; star/star: +1. Canonical20 scores unchanged.\n");
		out.append("# Statistical E-value parameters are not provided.\n");
		out.append("# BBTools import validates extras but retains its existing scaled-BLOSUM X and invalid-stop runtime policy.\n ");
		for(int a=0; a<OUTPUT_SYMBOLS.length(); a++){out.append(' ').append(OUTPUT_SYMBOLS.charAt(a));}
		out.nl();
		int minimum=Integer.MAX_VALUE;
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){minimum=Math.min(minimum, matrix.score((byte)a, (byte)b));}
		}
		for(int a=0; a<OUTPUT_SYMBOLS.length(); a++){
			out.append(OUTPUT_SYMBOLS.charAt(a));
			for(int b=0; b<OUTPUT_SYMBOLS.length(); b++){
				final int score;
				if(a==23 || b==23){score=a==23 && b==23 ? 1 : minimum;}
				else if(a<20 && b<20){score=matrix.score((byte)code(a), (byte)code(b));}
				else{score=expected(matrix, weights[Math.min(a, b)], weights[Math.max(a, b)]);}
				out.append(' ').append(score);
			}
			out.nl();
		}
		return out;
	}

	/** Validates the probability distribution and builds canonical/B/Z/X weights.
	 * Zero individual frequencies are allowed; zero D+N or E+Q makes a requested
	 * ambiguity distribution undefined and is rejected rather than guessed.
	 */
	private static double[][] weights(double[] background){
		require(background!=null && background.length==20, "Background requires20 encoded residue probabilities");
		double sum=0;
		for(double p:background){
			require(Double.isFinite(p) && p>=0 && p<=1, "Background probabilities must be finite in [0,1]");
			sum+=p;
		}
		require(Math.abs(sum-1)<=1e-9, "Background probabilities must sum to one");
		final double[][] weights=new double[23][20];
		for(int a=0; a<20; a++){weights[a][code(a)]=1;}
		for(int a=0; a<20; a++){weights[22][a]=background[a]/sum;}
		mixture(weights[20], background, 'D', 'N');
		mixture(weights[21], background, 'E', 'Q');
		return weights;
	}

	/** Fills a two-residue conditional distribution after checking its denominator. */
	private static void mixture(double[] weights, double[] background, char first, char second){
		final int a=AminoAcid.acidToNumberExtended[first], b=AminoAcid.acidToNumberExtended[second];
		final double total=background[a]+background[b];
		require(total>0, "Undefined background distribution for "+first+"/"+second);
		weights[a]=background[a]/total;
		weights[b]=background[b]/total;
	}

	/** Rounds a symmetric weighted score once, using the producer's half-away rule. */
	private static int expected(ProteinSubstitutionMatrix matrix, double[] first, double[] second){
		double value=0;
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){value+=first[a]*second[b]*matrix.score((byte)a, (byte)b);}
		}
		assert(Double.isFinite(value) && Math.abs(value)<=101) : "Normalized mixtures preserve bounded matrix scores";
		return value<0 ? (int)(value-.5) : (int)(value+.5);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Conversion CLI        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads a matrix and optional authoritative normalized monomer state before
	 * opening output. Conversion never writes an intermediate file or refits scores.
	 */
	private static void convert(String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		String format="ncbi", name=null, background=null, band=null;
		int scale=0, gap=0;
		for(String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase(Locale.ROOT), value=split.length>1 ? split[1] : null;
			if(key.equals("outformat")){format=value;}
			else if(key.equals("name")){name=value;}
			else if(key.equals("scale")){scale=Integer.parseInt(value);}
			else if(key.equals("gap")){gap=Integer.parseInt(value);}
			else if(key.equals("background")){background=value;}
			else if(key.equals("band")){band=value;}
			else{require(parser.parse(arg, key, value), "Unknown parameter: "+arg);}
		}
		require(parser.in1!=null && parser.out1!=null && !parser.append, "Require in/out; append is unsupported");
		require("native".equals(format) || "ncbi".equals(format), "outformat must be native or ncbi");
		require(Tools.testInputFiles(false, false, parser.in1, background), "Input matrix/background is unreadable");
		require(Tools.testOutputFiles(parser.overwrite, false, false, parser.out1), "Output exists or is unwritable");
		require(Tools.testForDuplicateFiles(false, parser.in1, background, parser.out1), "Inputs and output must differ");
		protectInput(parser.in1, parser.out1);
		if(background!=null){protectInput(background, parser.out1);}
		final ProteinSubstitutionMatrix matrix=read(parser.in1, name, scale, gap);
		final ByteBuilder output;
		if(format.equals("native")){
			require(background==null && band==null, "Native export does not consume background or band");
			output=nativeText(matrix);
		}else{
			require(background!=null && band!=null && !band.isEmpty(), "NCBI export requires background=PROBABILITIES and band=NAME");
			if(ProteinSubstitutionMatrix.POSSUM55_NAME.equals(matrix.name)){
				require(band.equals("25to55"), "POSSUM55 export requires the accepted25to55 background band");
			}
			final FamilySubstitutionCountBlend.Table table=FamilySubstitutionCountBlend.readProbabilities(
				Paths.get(background), false, new String[]{band}).get(band);
			require(table!=null && table.n==20, "Missing monomer background band");
			final double[] raw=table.endpointFrequencies(), frequencies=new double[20];
			final String[] labels=FamilySubstitutionCountBlend.labels(false);
			for(int i=0; i<20; i++){frequencies[AminoAcid.acidToNumberExtended[labels[i].charAt(0)]]=raw[i];}
			final boolean possum=ProteinSubstitutionMatrix.POSSUM55_NAME.equals(matrix.name);
			if(possum){verifyPossumScores(matrix, table);}
			output=ncbiText(matrix, frequencies, background+"; band="+band+"; (query+reference)/(2*total)", possum);
		}
		final ByteStreamWriter writer=new ByteStreamWriter(parser.out1, parser.overwrite, false, true);
		writer.start();
		try{writer.print(output);}
		finally{require(!writer.poisonAndWait(), "Matrix output failed: "+parser.out1);}
	}

	/** Checks all400 fitted scores before attaching the named derivation formula.
	 * The caller has selected25to55. This proves the supplied matrix/state pairing,
	 * not that arbitrary input data came from the project's accepted R5 population.
	 */
	private static void verifyPossumScores(ProteinSubstitutionMatrix matrix, FamilySubstitutionCountBlend.Table table){
		require(matrix.scale==4 && matrix.linearGap==16 && table.n==20,
			"POSSUM55 derivation requires20 residues, native scale4 and gap16");
		final double[] odds=table.odds();
		final String[] labels=FamilySubstitutionCountBlend.labels(false);
		for(int a=0; a<20; a++){
			for(int b=0; b<20; b++){
				final double value=8*Math.log(odds[a*20+b])/Math.log(2);
				require(Double.isFinite(value) && Math.abs(value)<=101, "Invalid POSSUM55 probability-derived score");
				final int rounded=value<0 ? (int)(value-.5) : (int)(value+.5);
				final byte ca=(byte)AminoAcid.acidToNumberExtended[labels[a].charAt(0)];
				final byte cb=(byte)AminoAcid.acidToNumberExtended[labels[b].charAt(0)];
				require(matrix.score(ca, cb)==rounded,
					"POSSUM55 matrix/background mismatch at "+labels[a]+"/"+labels[b]);
			}
		}
	}

	/** Checks actual existing-file identity as well as normalized paths before output can overwrite input. */
	private static void protectInput(String input, String output){
		final Path in=Paths.get(input).toAbsolutePath().normalize(), out=Paths.get(output).toAbsolutePath().normalize();
		require(!in.equals(out), "Matrix conversion cannot overwrite its input");
		try{require(!Files.exists(out) || !Files.isSameFile(in, out), "Output aliases an input file");}
		catch(java.io.IOException e){throw new IllegalArgumentException("Cannot verify input/output separation", e);}
	}

	/** Stops before accepting malformed matrix state or emitting misleading score units. */
	private static void require(boolean valid, String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	private static final String CANONICAL="ARNDCQEGHILKMFPSTWYV";
	private static final String OUTPUT_SYMBOLS=CANONICAL+"BZX*";
	private static final String INPUT_SYMBOLS=OUTPUT_SYMBOLS+"J";
}
