package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Derives symmetric, rounded log-odds matrices from pooled monomer counts.
 * The arithmetic order and binary64 mass tolerance match the accepted
 * derive_matrix() recipe; integer weights have already been applied upstream.
 *
 * @author Sayu
 */
public final class FamilySubstitutionMatrixDerive{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static command-line tool; instances carry no state. */
	private FamilySubstitutionMatrixDerive(){}

	/** Parses standard BBTools arguments and closes the preliminary output stream. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilySubstitutionMatrixDerive.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes one matrix and population record per requested band into a fresh directory. */
	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final Path input=Paths.get(required(options, "in")).toAbsolutePath().normalize();
		final Path output=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		final String name=options.containsKey("name") ? options.get("name") : ProteinSubstitutionMatrix.DEFAULT_NAME;
		require(name.indexOf('\t')<0 && name.indexOf('\n')<0 && name.indexOf('\r')<0,
			"Matrix name must fit one TSV field");
		require(Files.isRegularFile(input) && !Files.isSymbolicLink(input), "Missing regular input: "+input);
		require(!Files.exists(output) && Files.isDirectory(output.getParent()),
			"Output must be a fresh directory under an existing parent: "+output);
		final String format=options.containsKey("format") ? options.get("format") : "counts";
		require(format.equals("counts") || format.equals("normalized_probabilities"), "format must be counts or normalized_probabilities");
		if(format.equals("normalized_probabilities")){
			deriveProbabilities(input, output, name, FamilySubstitutionCountBlend.bands(required(options, "bands")));
			return;
		}
		final LinkedHashMap<String,Band> bands=new LinkedHashMap<String,Band>();
		for(String band:required(options, "bands").split(",", -1)){
			require(band.matches("[A-Za-z0-9][A-Za-z0-9_.-]*"), "Unsafe or empty band: "+band);
			require(bands.put(band, new Band(band))==null, "Duplicate requested band: "+band);
		}
		read(input, bands);
		for(Band band:bands.values()){band.finish();}
		Files.createDirectory(output);
		for(Band band:bands.values()){
			final ByteBuilder matrix=new ByteBuilder("#bbgap_substitution_matrix\t1\nname\t");
			matrix.append(name).append("\nscale\t4\ngap\t16\n");
			for(int i=0; i<SYMBOLS; i++){
				for(int j=0; j<SYMBOLS; j++){
					matrix.append("score").tab().append(ALPHABET.charAt(i)).tab()
						.append(ALPHABET.charAt(j)).tab().append(band.scores[i][j]).nl();
				}
			}
			write(output.resolve(band.name+".tsv"), matrix.toString());
			write(output.resolve(band.name+".population"),
				band.name+"\t"+band.alignmentText+"\t"+band.totalText+"\n");
		}
	}

	/** Derives matrices from recursive fractional state without pretending it contains alignment counts. */
	private static void deriveProbabilities(Path input, Path output, String name, String[] bands) throws java.io.IOException{
		final java.util.LinkedHashMap<String,FamilySubstitutionCountBlend.Table> tables=
			FamilySubstitutionCountBlend.readProbabilities(input, false, bands);
		Files.createDirectory(output);
		for(String band:bands){
			final FamilySubstitutionCountBlend.Table table=tables.get(band);
			final double[] odds=table.odds();
			final ByteBuilder matrix=new ByteBuilder("#bbgap_substitution_matrix\t1\nname\t");
			matrix.append(name).append("\nscale\t4\ngap\t16\n");
			for(int i=0; i<SYMBOLS; i++){
				for(int j=0; j<SYMBOLS; j++){
					require(odds[i*SYMBOLS+j]>0, "Zero observed normalized cell in band "+band);
					final double value=8*Math.log(odds[i*SYMBOLS+j])/LOG2;
					require(Double.isFinite(value) && value>Integer.MIN_VALUE+1.0 && value<Integer.MAX_VALUE-1.0, "Invalid normalized matrix score");
					final int score=value<0 ? (int)(value-.5) : (int)(value+.5);
					matrix.append("score").tab().append(ALPHABET.charAt(i)).tab().append(ALPHABET.charAt(j)).tab().append(score).nl();
				}
			}
			write(output.resolve(band+".tsv"), matrix.toString());
			write(output.resolve(band+".population"), "#bbgap_substitution_probability_population\t1\nband\tprobability_mass\n"
				+band+"\t"+Double.toString(table.total())+"\n");
		}
	}

	/** Reads the selected ordered 20-by-20 tables, retaining repeated endpoint masses
	 * and population text while rejecting duplicate cells or inconsistent denominators. */
	private static void read(Path input, LinkedHashMap<String,Band> bands){
		final ByteFile reader=ByteFile.makeByteFile(input.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			final byte[] first=reader.nextLine();
			require(first!=null && HEADER.equals(new String(first, StandardCharsets.US_ASCII)), "Unexpected pooled monomer header");
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				final Band band=bands.get(parser.parseString(0));
				if(band==null){continue;}
				require(parser.terms()==8, "Selected monomer row requires eight fields: "+band.name);
				final int i=residue(parser.parseString(1)), j=residue(parser.parseString(2));
				require(!band.seen[i][j], "Duplicate cell in band "+band.name+": "+i+","+j);
				band.seen[i][j]=true;
				band.counts[i][j]=count(parser.parseString(3));
				band.query[i][j]=count(parser.parseString(4));
				band.reference[i][j]=count(parser.parseString(5));
				final String alignmentText=parser.parseString(6), totalText=parser.parseString(7);
				final double alignments=count(alignmentText), total=count(totalText);
				if(band.cells==0){
					band.alignments=alignments;
					band.total=total;
					band.alignmentText=alignmentText;
					band.totalText=totalText;
				}else{
					require(band.alignments==alignments && band.total==total, "Mismatched denominators in band "+band.name);
				}
				band.cells++;
			}
		}finally{require(!reader.close(), "Failed reading monomer counts: "+input);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One selected band with directed counts, endpoint witnesses and rounded scores. */
	private static final class Band{
		/** Allocates the fixed-size tables for one output band. */
		Band(String name_){name=name_;}

		/** Checks complete-cell and mass invariants, then reproduces the recipe's
		 * binary64 symmetric log-odds calculation and half-away integer rounding. */
		void finish(){
			require(cells==SYMBOLS*SYMBOLS && alignments>=1 && total>=1,
				"Band needs 400 cells and positive denominators: "+name);
			// Deliberately accumulate as binary64 in the recipe's residue order.
			// Exact long summation would change the accepted large weighted-count arithmetic.
			final double[] endpoints=new double[SYMBOLS];
			double sum=0;
			for(int i=0; i<SYMBOLS; i++){
				double q=0, r=0;
				for(int j=0; j<SYMBOLS; j++){
					q+=counts[i][j];
					r+=counts[j][i];
				}
				for(int j=0; j<SYMBOLS; j++){
					require(nearly(query[i][j], q) && nearly(reference[j][i], r),
						"Mismatched endpoint mass in band "+name);
				}
				endpoints[i]=q+r;
				sum+=endpoints[i];
			}
			require(nearly(sum, 2*total), "Endpoint sum differs from twice the total in band "+name);
			for(int i=0; i<SYMBOLS; i++){
				for(int j=0; j<SYMBOLS; j++){
					final double observed=i==j ? counts[i][i] : counts[i][j]+counts[j][i];
					require(observed>0, "Zero observed cell in band "+name+": "+i+","+j);
					final double expected=(endpoints[i]/(2*total))*(endpoints[j]/(2*total))*(i==j ? 1 : 2);
					final double value=8*Math.log((observed/total)/expected)/LOG2;
					require(Double.isFinite(value) && value>Integer.MIN_VALUE+1.0 && value<Integer.MAX_VALUE-1.0,
						"Invalid log-odds score in band "+name);
					scores[i][j]=value<0 ? (int)(value-.5) : (int)(value+.5);
				}
			}
		}

		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/

		final String name;
		final double[][] counts=new double[SYMBOLS][SYMBOLS];
		final double[][] query=new double[SYMBOLS][SYMBOLS];
		final double[][] reference=new double[SYMBOLS][SYMBOLS];
		final boolean[][] seen=new boolean[SYMBOLS][SYMBOLS];
		final int[][] scores=new int[SYMBOLS][SYMBOLS];
		int cells;
		double alignments, total;
		String alignmentText, totalText;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Matches the accepted weighted-count recipe's absolute endpoint tolerance. */
	private static boolean nearly(double a, double b){return Math.abs(a-b)<=4096;}

	/** Resolves one standard amino acid in the fixed matrix output order. */
	private static int residue(String text){
		final int index=text.length()==1 ? ALPHABET.indexOf(text.charAt(0)) : -1;
		require(index>=0, "Bad amino acid: "+text);
		return index;
	}

	/** Validates unsigned decimal count text and converts directly to recipe binary64. */
	private static double count(String text){
		require(text.length()>0, "Empty count");
		for(int i=0; i<text.length(); i++){
			final char c=text.charAt(i);
			require(c>='0' && c<='9', "Nonnumeric count: "+text);
		}
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value), "Nonfinite count: "+text);
		return value;
	}

	/** Reads supported flag=value options, rejecting unknown or repeated keys. */
	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> out=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && pair[1].length()>0, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(key.equals("in") || key.equals("out") || key.equals("bands") || key.equals("name") || key.equals("format"), "Unknown option: "+key);
			require(out.put(key, pair[1])==null, "Duplicate option: "+key);
		}
		return out;
	}

	/** Returns a required nonempty option or explains the missing argument. */
	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		require(value!=null && value.length()>0, "Missing option: "+key);
		return value;
	}

	/** Creates one ASCII matrix or population file without overwriting existing bytes. */
	private static void write(Path path, String text) throws java.io.IOException{
		Files.write(path, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Stops on invalid input or a violated scientific invariant. */
	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final String ALPHABET="ARNDCEQGHILKMFPSTWYV";
	private static final int SYMBOLS=20;
	private static final double LOG2=Math.log(2);
	private static final String HEADER="band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs";
}
