package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Normalizes and blends directed monomer/dimer counts before odds estimation.
 * Versioned probability output is also the previous-input state for the next round.
 * @author Sayu
 */
public final class FamilySubstitutionCountBlend{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Static tool and shared probability-table API. */
	private FamilySubstitutionCountBlend(){}

	/** Runs the explicit weighted blend and releases the standard argument stream. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, FamilySubstitutionCountBlend.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Blends both kinds independently at each band; weights are relative, never count multipliers. */
	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final double newWeight=weight(required(options, "newweight"));
		final double previousWeight=weight(required(options, "previousweight"));
		require(Double.isFinite(newWeight+previousWeight) && newWeight+previousWeight>0, "Weights need a finite positive sum");
		final String[] bands=bands(required(options, "bands"));
		final Path output=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		require(!Files.exists(output) && Files.isDirectory(output.getParent()), "Output must be a fresh directory: "+output);
		final LinkedHashMap<String,Table> monomers=blendInputs(options, false, bands, newWeight, previousWeight);
		final LinkedHashMap<String,Table> dimers=blendInputs(options, true, bands, newWeight, previousWeight);
		Files.createDirectory(output);
		Files.createDirectory(output.resolve("monomer"));
		Files.createDirectory(output.resolve("dimer"));
		write(output.resolve("monomer/probabilities.tsv"), false, monomers);
		write(output.resolve("dimer/probabilities.tsv"), true, dimers);
		final ByteBuilder manifest=new ByteBuilder("#bbgap_substitution_count_blend\t1\n");
		for(String key:new String[]{"newmonomers", "newdimers", "previousmonomers", "previousdimers", "bands", "newweight", "previousweight"}){
			manifest.append(key).tab().append(options.get(key)).nl();
		}
		manifest.append("new_coefficient\t").append(Double.toString(newWeight/(newWeight+previousWeight))).nl()
			.append("previous_coefficient\t").append(Double.toString(previousWeight/(newWeight+previousWeight))).nl()
			.append("normalization\tindependent unit total per input, band and monomer/dimer kind\n")
			.append("previous_state\tserialized previous blend, once blending has begun\n");
		Files.write(output.resolve("manifest.tsv"), manifest.toString().getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE_NEW);
		Files.write(output.resolve("PASS"), "FAMILY_SUBSTITUTION_COUNT_BLEND_PASS\n".getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Loads the two independently normalized input populations for one alphabet size. */
	private static LinkedHashMap<String,Table> blendInputs(HashMap<String,String> options, boolean dimer,
			String[] bands, double newWeight, double previousWeight){
		final String suffix=dimer ? "dimers" : "monomers";
		final LinkedHashMap<String,Table> current=read(Paths.get(required(options, "new"+suffix)), dimer, bands);
		final LinkedHashMap<String,Table> previous=read(Paths.get(required(options, "previous"+suffix)), dimer, bands);
		final LinkedHashMap<String,Table> result=new LinkedHashMap<String,Table>();
		for(String band:bands){result.put(band, blend(current.get(band), previous.get(band), newWeight, previousWeight));}
		return result;
	}

	/** Applies relative coefficients to separately unit-normalized tables, including a prior blend. */
	public static Table blend(Table current, Table previous, double newWeight, double previousWeight){
		require(current.n==previous.n && current.band.equals(previous.band), "Blend populations must have the same band and symbols");
		require(Double.isFinite(newWeight) && Double.isFinite(previousWeight) && newWeight>=0 && previousWeight>=0
			&& Double.isFinite(newWeight+previousWeight) && newWeight+previousWeight>0, "Invalid blend weights");
		final double newCoefficient=newWeight/(newWeight+previousWeight);
		final double previousCoefficient=previousWeight/(newWeight+previousWeight);
		final Table out=new Table(current.band, current.n);
		for(int i=0; i<out.values.length; i++){
			out.values[i]=previousCoefficient*(previous.values[i]/previous.total())+newCoefficient*(current.values[i]/current.total());
		}
		out.finish();
		require(Math.abs(out.total()-1)<=TOLERANCE, "Blended probability mass must remain one");
		return out;
	}

	/** Reads exact raw pooled counts or the explicitly versioned recursive state, never an old 50:50 schema. */
	public static LinkedHashMap<String,Table> read(Path input, boolean dimer, String[] selectedBands){
		return read(input, dimer, selectedBands, false);
	}

	/** Requires versioned normalized state for scientific consumers selecting the probability route. */
	public static LinkedHashMap<String,Table> readProbabilities(Path input, boolean dimer, String[] selectedBands){
		return read(input, dimer, selectedBands, true);
	}

	/** Shares schema/mass validation while keeping the explicit probability route fail-closed. */
	private static LinkedHashMap<String,Table> read(Path input, boolean dimer, String[] selectedBands, boolean requireProbabilities){
		require(Files.isRegularFile(input) && !Files.isSymbolicLink(input), "Missing regular count/state input: "+input);
		final String[] labels=labels(dimer);
		final LinkedHashMap<String,ReadState> states=new LinkedHashMap<String,ReadState>();
		for(String band:selectedBands){require(states.put(band, new ReadState(new Table(band, labels.length)))==null, "Duplicate selected band");}
		final ByteFile reader=ByteFile.makeByteFile(input.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		boolean probabilities=false;
		try{
			final String first=text(reader.nextLine());
			probabilities=MAGIC.equals(first);
			require(!requireProbabilities || probabilities, "Explicit normalized_probabilities format requires versioned probability state");
			if(probabilities){
				require(("kind\t"+(dimer ? "dimer" : "monomer")).equals(text(reader.nextLine())), "Probability-state kind mismatch");
				require(STATE_HEADER.equals(text(reader.nextLine())), "Unexpected probability-state columns");
			}else{require((dimer ? DIMER_HEADER : MONOMER_HEADER).equals(first), "Unexpected raw pooled count schema");}
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				final ReadState state=states.get(parser.parseString(0));
				if(state==null){continue;}
				require(parser.terms()==(probabilities ? 6 : dimer ? 4 : 8), "Incorrect count/state field count");
				final int cell=state.cells++, n=labels.length, i=cell/n, j=cell%n;
				require(cell<n*n && parser.termEquals(labels[i], 1) && parser.termEquals(labels[j], 2),
					"Every selected band must contain the complete canonical directed table exactly once");
				if(probabilities){
					state.table.values[cell]=probability(parser.parseString(3));
					state.marginal(state.query, i, probability(parser.parseString(4)));
					state.marginal(state.reference, j, probability(parser.parseString(5)));
				}else{
					final long count=parser.parseLong(3);
					require(count>=0, "Raw counts must be nonnegative");
					state.table.values[cell]=count;
					state.rawTotal=Math.addExact(state.rawTotal, count);
					state.rawQuery[i]=Math.addExact(state.rawQuery[i], count);
					state.rawReference[j]=Math.addExact(state.rawReference[j], count);
					if(!dimer){state.monomerWitness(parser, i, j);}
				}
			}
		}finally{require(!reader.close(), "Failed reading count/state input: "+input);}
		final LinkedHashMap<String,Table> result=new LinkedHashMap<String,Table>();
		for(ReadState state:states.values()){
			require(state.cells==labels.length*labels.length, "Incomplete selected band: "+state.table.band);
			state.table.finish();
			if(probabilities){
				require(Math.abs(state.table.total()-1)<=TOLERANCE, "Saved probability mass must be one");
				for(int i=0; i<labels.length; i++){
					require(Math.abs(state.query[i]-state.table.queryMarginals[i])<=TOLERANCE
						&& Math.abs(state.reference[i]-state.table.referenceMarginals[i])<=TOLERANCE, "Saved marginal disagrees with directed probabilities");
				}
			}else{
				require(state.rawTotal>0, "Raw population must have positive mass");
				if(!dimer){
					require(state.rawTotal==state.declaredTotal, "Raw monomer total differs from cell sum");
					for(int i=0; i<labels.length; i++){
						require(state.rawQuery[i]==state.queryWitness[i] && state.rawReference[i]==state.referenceWitness[i], "Raw endpoint mass mismatch");
					}
				}
				// Preserve the exact integer total before conversion, rather than summing rounded large doubles.
				state.table.total=state.rawTotal;
			}
			result.put(state.table.band, state.table);
		}
		return result;
	}

	/** Writes losslessly round-trippable probabilities and marginals for recursive reuse. */
	public static void write(Path path, boolean dimer, LinkedHashMap<String,Table> tables){
		require(!Files.exists(path), "Refusing to overwrite probability state: "+path);
		final String[] labels=labels(dimer);
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		final ByteBuilder row=new ByteBuilder();
		writer.start();
		try{
			writer.println(MAGIC);
			writer.println("kind\t"+(dimer ? "dimer" : "monomer"));
			writer.println(STATE_HEADER);
			for(Table table:tables.values()){
				require(table.n==labels.length && Math.abs(table.total()-1)<=TOLERANCE, "Only normalized probability state may be serialized");
				for(int i=0; i<table.n; i++){
					for(int j=0; j<table.n; j++){
						row.clear().append(table.band).tab().append(labels[i]).tab().append(labels[j]).tab()
							.append(Double.toString(table.values[i*table.n+j])).tab()
							.append(Double.toString(table.queryMarginals[i])).tab().append(Double.toString(table.referenceMarginals[j]));
						writer.println(row);
					}
				}
			}
		}finally{require(!writer.poisonAndWait(), "Probability-state write failed: "+path);}
	}

	/** Returns the established monomer order or the matching canonical 400-dimer dictionary. */
	public static String[] labels(boolean dimer){
		if(dimer){return FamilyDimerSubstitutionTable.dimerLabels();}
		final String[] out=new String[20];
		for(int i=0; i<out.length; i++){out[i]=ALPHABET.substring(i, i+1);}
		return out;
	}

	/** Validates a nonempty band list without filesystem-derived keys. */
	public static String[] bands(String csv){
		final String[] out=csv.split(",", -1);
		for(String band:out){require(band.matches("[A-Za-z0-9][A-Za-z0-9_.-]*"), "Invalid band: "+band);}
		return out;
	}

	/** Converts one finite, nonnegative state value; normalization is checked after the table. */
	private static double probability(String text){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>=0 && value<=1+TOLERANCE, "Invalid probability/marginal: "+text);
		return value;
	}

	/** Validates a finite nonnegative relative blend weight. */
	private static double weight(String text){
		final double value=Double.parseDouble(text);
		require(Double.isFinite(value) && value>=0, "Invalid blend weight: "+text);
		return value;
	}

	/** Reads an ASCII control/header line without accepting missing fields. */
	private static String text(byte[] line){require(line!=null, "Missing count/state header");return new String(line, StandardCharsets.US_ASCII);}

	/** Parses only the explicit blend parameters, rejecting duplicate keys. */
	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> out=new HashMap<String,String>();
		for(String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && pair[1].length()>0, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(Arrays.asList("newmonomers", "newdimers", "previousmonomers", "previousdimers", "bands", "newweight", "previousweight", "out").contains(key)
				&& out.put(key, pair[1])==null, "Unknown or repeated option: "+key);
		}
		return out;
	}

	/** Returns a required nonempty option. */
	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		require(value!=null && value.length()>0, "Missing option: "+key);
		return value;
	}

	/** Stops before emitting incorrect scientific state. */
	private static void require(boolean condition, String message){if(!condition){throw new IllegalArgumentException(message);}}

	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Directed masses and double marginals; raw integer totals are retained until normalization. */
	public static final class Table{
		/** Allocates one complete directed table. */
		private Table(String band_, int n_){
			band=band_;
			n=n_;
			values=new double[n*n];
			queryMarginals=new double[n];
			referenceMarginals=new double[n];
		}

		/** Accumulates masses in canonical row-major order. */
		private void finish(){
			Arrays.fill(queryMarginals, 0);Arrays.fill(referenceMarginals, 0);total=0;
			for(int i=0; i<n; i++){
				for(int j=0; j<n; j++){
					final double value=values[i*n+j];
					require(Double.isFinite(value) && value>=0, "Invalid directed mass");
					queryMarginals[i]+=value;referenceMarginals[j]+=value;total+=value;
				}
			}
			require(Double.isFinite(total) && total>0, "Directed table needs finite positive total mass");
		}

		/** Returns the actual input mass used for independent normalization. */
		public double total(){return total;}

		/** Returns symmetric endpoint probabilities without integer-rescaling fractional counts. */
		public double[] endpointFrequencies(){
			final double[] out=new double[n];
			for(int i=0; i<n; i++){out[i]=(queryMarginals[i]+referenceMarginals[i])/(2*total);}
			return out;
		}

		/** Computes symmetric observed/expected odds before any matrix rounding or alphabet reduction. */
		public double[] odds(){
			final double[] frequency=endpointFrequencies(), out=new double[n*n];
			for(int i=0; i<n; i++){
				for(int j=i; j<n; j++){
					final double observed=(i==j ? values[i*n+i] : values[i*n+j]+values[j*n+i])/total;
					final double expected=(frequency[i]*frequency[j])*(i==j ? 1 : 2);
					final double value=expected>0 ? observed/expected : 0;
					require(Double.isFinite(value) && value>=0, "Invalid normalized odds");
					out[i*n+j]=out[j*n+i]=value;
				}
			}
			return out;
		}

		final String band;
		public final int n;
		public final double[] values, queryMarginals, referenceMarginals;
		private double total;
	}

	/** Validation witnesses for raw counts or saved probability marginals. */
	private static final class ReadState{
		/** Allocates only per-symbol marginal witnesses, not another directed table. */
		ReadState(Table table_){
			table=table_;
			final int n=table.n;
			query=new double[n];
			reference=new double[n];
			Arrays.fill(query, -1);
			Arrays.fill(reference, -1);
			rawQuery=new long[n];
			rawReference=new long[n];
			queryWitness=new long[n];
			referenceWitness=new long[n];
			Arrays.fill(queryWitness, -1);
			Arrays.fill(referenceWitness, -1);
		}

		/** Requires repeated serialized marginals to agree exactly before mass checks. */
		void marginal(double[] witnesses, int index, double value){
			if(witnesses[index]<0){witnesses[index]=value;}
			else{require(Double.doubleToLongBits(witnesses[index])==Double.doubleToLongBits(value), "Inconsistent repeated marginal");}
		}

		/** Retains and checks the exact raw monomer endpoint/population witnesses. */
		void monomerWitness(LineParser1 parser, int i, int j){
			final long q=parser.parseLong(4), r=parser.parseLong(5), a=parser.parseLong(6), t=parser.parseLong(7);
			require(q>=0 && r>=0 && a>0 && t>0, "Invalid monomer population witnesses");
			if(queryWitness[i]<0){queryWitness[i]=q;}else{require(queryWitness[i]==q, "Inconsistent query marginal");}
			if(referenceWitness[j]<0){referenceWitness[j]=r;}else{require(referenceWitness[j]==r, "Inconsistent reference marginal");}
			if(cells==1){alignments=a;declaredTotal=t;}
			else{require(alignments==a && declaredTotal==t, "Inconsistent raw denominators");}
		}

		final Table table;
		final double[] query, reference;
		final long[] rawQuery, rawReference, queryWitness, referenceWitness;
		long rawTotal, declaredTotal, alignments;
		int cells;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	public static final String MAGIC="#bbgap_substitution_probabilities\t1";
	private static final String STATE_HEADER="band\tquery_symbol\tconsensus_symbol\tprobability\tquery_marginal\tconsensus_marginal";
	private static final String MONOMER_HEADER="band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs";
	private static final String DIMER_HEADER="band\tquery_dimer\tconsensus_dimer\tcount";
	private static final String ALPHABET="ARNDCEQGHILKMFPSTWYV";
	private static final double TOLERANCE=1e-10;
}
