package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Pools independently counted SHORT/LONG family substitution tables with integer
 * population weights before any matrix or alphabet normalization.
 *
 * @author Yelan
 */
public final class FamilyWeightedSubstitutionTables{

	private FamilyWeightedSubstitutionTables(){}

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, FamilyWeightedSubstitutionTables.class, false);
			try{
				if(pp.args.length==1 && pp.args[0].equalsIgnoreCase("test")){test();}
				else{run(pp.args);}
			}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	public static void run(String[] args) throws Exception{
		final HashMap<String,String> options=parse(args);
		final Path[] monomers=paths(required(options, "monomers"));
		final Path[] dimers=paths(required(options, "dimers"));
		final long[] weights=weights(required(options, "weights"));
		final Path out=Paths.get(required(options, "out")).toAbsolutePath().normalize();
		require(monomers.length==dimers.length && monomers.length==weights.length && weights.length>0,
			"monomers, dimers and weights must describe the same nonempty source list");
		require(Files.isDirectory(out.getParent()) && !Files.exists(out), "Output directory must be fresh");

		MonomerTable monomer=null;
		DimerTable dimer=null;
		for(int i=0; i<weights.length; i++){
			monomer=MonomerTable.add(monomer, readMonomer(monomers[i]), weights[i]);
			dimer=DimerTable.add(dimer, readDimer(dimers[i]), weights[i]);
		}
		require(Arrays.equals(monomer.bands, dimer.bands), "Monomer and dimer band labels differ");

		Files.createDirectory(out);
		Files.createDirectory(out.resolve("monomer"));
		Files.createDirectory(out.resolve("dimer"));
		writeMonomer(out.resolve("monomer/ordered_counts.tsv"), monomer);
		writeDimer(out.resolve("dimer"), dimer);
		writeManifest(out.resolve("manifest.tsv"), monomers, dimers, weights);
		write(out.resolve("PASS"), "FAMILY_WEIGHTED_SUBSTITUTION_TABLES_PASS sources="+weights.length+
			" bands="+monomer.bands.length+"\n");
	}

	private static MonomerTable readMonomer(Path path){
		final LinkedHashMap<String,long[][]> counts=new LinkedHashMap<String,long[][]>();
		final HashMap<String,long[]> alignments=new HashMap<String,long[]>();
		final HashMap<String,boolean[][]> seen=new HashMap<String,boolean[][]>();
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile(path.toString(), false);
		try{
			require(equals(reader.nextLine(), MONOMER_HEADER), "Unexpected monomer schema: "+path);
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				require(parser.terms()==8, "Monomer rows require eight columns: "+path);
				final String band=parser.parseString(0);
				final int q=residue(parser.parseString(1)), r=residue(parser.parseString(2));
				long[][] table=counts.get(band);
				if(table==null){
					table=new long[20][20];
					counts.put(band, table);
					seen.put(band, new boolean[20][20]);
					alignments.put(band, new long[]{Long.parseLong(parser.parseString(6)), Long.parseLong(parser.parseString(7))});
				}
				require(!seen.get(band)[q][r], "Duplicate monomer cell: "+band+" "+q+","+r);
				seen.get(band)[q][r]=true;
				table[q][r]=parseCount(parser.parseString(3));
				require(alignments.get(band)[0]==Long.parseLong(parser.parseString(6)) &&
					alignments.get(band)[1]==Long.parseLong(parser.parseString(7)), "Repeated monomer denominators differ: "+band);
			}
		}finally{require(!reader.close(), "Failed reading monomer table: "+path);}
		final MonomerTable out=new MonomerTable(counts.keySet().toArray(new String[0]));
		for(int b=0; b<out.bands.length; b++){
			final String band=out.bands[b];
			final boolean[][] used=seen.get(band);
			long total=0;
			for(int q=0; q<20; q++){
				for(int r=0; r<20; r++){
					require(used[q][r], "Missing monomer cell: "+band+" "+q+","+r);
					out.counts[b][q][r]=counts.get(band)[q][r];
					total=Math.addExact(total, out.counts[b][q][r]);
				}
			}
			out.alignments[b]=alignments.get(band)[0];
			require(total==alignments.get(band)[1], "Monomer aligned-pair denominator differs from cell sum: "+band);
		}
		return out;
	}

	private static DimerTable readDimer(Path path){
		final LinkedHashMap<String,long[]> counts=new LinkedHashMap<String,long[]>();
		final HashMap<String,boolean[]> seen=new HashMap<String,boolean[]>();
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile reader=ByteFile.makeByteFile(path.toString(), false);
		try{
			require(equals(reader.nextLine(), DIMER_HEADER), "Unexpected dimer schema: "+path);
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				parser.set(line);
				require(parser.terms()==4, "Dimer rows require four columns: "+path);
				final String band=parser.parseString(0);
				long[] table=counts.get(band);
				if(table==null){
					table=new long[DIMERS*DIMERS];
					counts.put(band, table);
					seen.put(band, new boolean[DIMERS*DIMERS]);
				}
				final int q=dimer(parser.parseString(1)), r=dimer(parser.parseString(2));
				final int cell=q*DIMERS+r;
				require(!seen.get(band)[cell], "Duplicate dimer cell: "+band+" "+q+","+r);
				seen.get(band)[cell]=true;
				table[cell]=parseCount(parser.parseString(3));
			}
		}finally{require(!reader.close(), "Failed reading dimer table: "+path);}
		final DimerTable out=new DimerTable(counts.keySet().toArray(new String[0]));
		for(int b=0; b<out.bands.length; b++){
			final String band=out.bands[b];
			for(int i=0; i<out.ordered[b].length; i++){
				require(seen.get(band)[i], "Missing dimer cell: "+band+" "+i);
				out.ordered[b][i]=counts.get(band)[i];
			}
		}
		return out;
	}

	private static void writeMonomer(Path path, MonomerTable table){
		final ByteStreamWriter writer=writer(path);
		final ByteBuilder row=new ByteBuilder();
		try{
			writer.println(MONOMER_HEADER_TEXT);
			for(int b=0; b<table.bands.length; b++){
				final long[] query=new long[20], consensus=new long[20];
				long total=0;
				for(int q=0; q<20; q++){for(int r=0; r<20; r++){
					query[q]=Math.addExact(query[q], table.counts[b][q][r]);
					consensus[r]=Math.addExact(consensus[r], table.counts[b][q][r]);
					total=Math.addExact(total, table.counts[b][q][r]);
				}}
				for(int q=0; q<20; q++){for(int r=0; r<20; r++){
					row.clear().append(table.bands[b]).tab().append((char)AminoAcid.numberToAcid[q]).tab()
						.append((char)AminoAcid.numberToAcid[r]).tab().append(table.counts[b][q][r]).tab()
						.append(query[q]).tab().append(consensus[r]).tab().append(table.alignments[b]).tab().append(total);
					writer.println(row);
				}}
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing monomer table: "+path);}
	}

	private static void writeDimer(Path root, DimerTable table){
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		ByteStreamWriter writer=writer(root.resolve("dictionary.tsv"));
		ByteBuilder row=new ByteBuilder();
		try{
			writer.println("dimer_id\tdimer");
			for(int i=0; i<labels.length; i++){
				row.clear().append(i).tab().append(labels[i]);
				writer.println(row);
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing dimer dictionary");}

		writer=writer(root.resolve("ordered_counts.tsv"));
		try{
			writer.println(DIMER_HEADER_TEXT);
			for(int b=0; b<table.bands.length; b++){
				for(int q=0; q<DIMERS; q++){for(int r=0; r<DIMERS; r++){
					row.clear().append(table.bands[b]).tab().append(labels[q]).tab().append(labels[r]).tab()
						.append(table.ordered[b][q*DIMERS+r]);
					writer.println(row);
				}}
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing ordered dimer table");}

		writer=writer(root.resolve("symmetric_tables.tsv"));
		try{
			writer.println("band\tdimer_i\tdimer_j\tcount\tendpoint_i\tendpoint_j\tobserved_frequency\texpected_frequency\todds");
			for(int b=0; b<table.bands.length; b++){
				final FamilyDimerSubstitutionTable.Table odds=new FamilyDimerSubstitutionTable.Table(table.ordered[b], DIMERS);
				require(odds.total>0, "Every pooled dimer band needs observations: "+table.bands[b]);
				for(int i=0; i<DIMERS; i++){for(int j=i; j<DIMERS; j++){
					final double pi=odds.endpoints[i]/(2.0*odds.total);
					final double pj=odds.endpoints[j]/(2.0*odds.total);
					row.clear().append(table.bands[b]).tab().append(labels[i]).tab().append(labels[j]).tab()
						.append(odds.counts[i*DIMERS+j]).tab().append(odds.endpoints[i]).tab().append(odds.endpoints[j]).tab()
						.append(Double.toString(odds.counts[i*DIMERS+j]/(double)odds.total)).tab()
						.append(Double.toString(i==j ? pi*pi : 2*pi*pj)).tab().append(Double.toString(odds.odds[i*DIMERS+j]));
					writer.println(row);
				}}
			}
		}finally{require(!writer.poisonAndWait(), "Failed writing symmetric dimer table");}
	}

	private static void writeManifest(Path path, Path[] monomers, Path[] dimers, long[] weights){
		final ByteBuilder out=new ByteBuilder("schema\tfamily_weighted_substitution_tables_v1\n");
		for(int i=0; i<weights.length; i++){
			out.append("source\t").append(i).tab().append("weight").tab().append(weights[i]).nl()
				.append("source\t").append(i).tab().append("monomers").tab().append(monomers[i].toString()).nl()
				.append("source\t").append(i).tab().append("dimers").tab().append(dimers[i].toString()).nl();
		}
		write(path, out.toString());
	}

	private static void test() throws Exception{
		final Path root=Files.createTempDirectory("family-weighted-substitutions-");
		final Path monomer1=root.resolve("m1.tsv"), monomer2=root.resolve("m2.tsv");
		final Path dimer1=root.resolve("d1.tsv"), dimer2=root.resolve("d2.tsv");
		write(monomer1, monomerFixture(1));
		write(monomer2, monomerFixture(2));
		write(dimer1, dimerFixture(1));
		write(dimer2, dimerFixture(2));
		run(new String[]{"monomers="+monomer1+","+monomer2, "dimers="+dimer1+","+dimer2,
			"weights=3,5", "out="+root.resolve("out")});
		final String monomer=new String(Files.readAllBytes(root.resolve("out/monomer/ordered_counts.tsv")), StandardCharsets.US_ASCII);
		final String dimer=new String(Files.readAllBytes(root.resolve("out/dimer/ordered_counts.tsv")), StandardCharsets.US_ASCII);
		require(monomer.contains("band0\tA\tA\t13\t"), "Monomer counts must be weighted before emission");
		require(dimer.contains("band0\tAA\tAA\t13\n"), "Dimer counts must be weighted before emission");
		require(Files.isRegularFile(root.resolve("out/dimer/dictionary.tsv")), "Dimer dictionary must be emitted");
		require(Files.isRegularFile(root.resolve("out/dimer/symmetric_tables.tsv")), "Dimer symmetric table must be emitted");
		System.out.println("FAMILY_WEIGHTED_SUBSTITUTION_TABLES_TEST_PASS");
	}

	private static String monomerFixture(long multiplier){
		final StringBuilder out=new StringBuilder(MONOMER_HEADER_TEXT).append('\n');
		final String alphabet="ARNDCEQGHILKMFPSTWYV";
		for(String band : new String[]{"band0", "band1"}){
			final long[][] counts=new long[20][20];
			final long[] query=new long[20], consensus=new long[20];
			long total=0;
			for(int q=0; q<20; q++){for(int r=0; r<20; r++){
				counts[q][r]=(q==r ? 1 : 2)*multiplier;
				query[q]+=counts[q][r];
				consensus[r]+=counts[q][r];
				total+=counts[q][r];
			}}
			for(int q=0; q<20; q++){for(int r=0; r<20; r++){
				out.append(band).append('\t').append(alphabet.charAt(q)).append('\t').append(alphabet.charAt(r)).append('\t')
					.append(counts[q][r]).append('\t').append(query[q]).append('\t').append(consensus[r]).append("\t7\t")
					.append(total).append('\n');
			}}
		}
		return out.toString();
	}

	private static String dimerFixture(long multiplier){
		final StringBuilder out=new StringBuilder(DIMER_HEADER_TEXT).append('\n');
		final String[] labels=FamilyDimerSubstitutionTable.dimerLabels();
		for(String band : new String[]{"band0", "band1"}){
			for(int q=0; q<DIMERS; q++){for(int r=0; r<DIMERS; r++){
				out.append(band).append('\t').append(labels[q]).append('\t').append(labels[r]).append('\t')
					.append((q==r ? 1 : 2)*multiplier).append('\n');
			}}
		}
		return out.toString();
	}

	private static HashMap<String,String> parse(String[] args){
		final HashMap<String,String> options=new HashMap<String,String>();
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			require(pair.length==2 && pair[1].length()>0, "Expected flag=value: "+arg);
			final String key=pair[0].toLowerCase(Locale.ROOT);
			require(OPTIONS.containsKey(key) && options.put(key, pair[1])==null, "Unknown or duplicate option: "+key);
		}
		return options;
	}

	private static Path[] paths(String csv){
		final String[] split=csv.split(",", -1);
		final Path[] out=new Path[split.length];
		for(int i=0; i<split.length; i++){
			require(split[i].length()>0, "Path lists cannot contain empty entries");
			out[i]=Paths.get(split[i]).toAbsolutePath().normalize();
			require(Files.isRegularFile(out[i]), "Missing input file: "+out[i]);
		}
		return out;
	}

	private static long[] weights(String csv){
		final String[] split=csv.split(",", -1);
		final long[] out=new long[split.length];
		for(int i=0; i<split.length; i++){
			out[i]=Long.parseLong(split[i]);
			require(out[i]>0, "Weights must be positive integers");
		}
		return out;
	}

	private static long parseCount(String s){
		final long value=Long.parseLong(s);
		require(value>=0, "Counts must be nonnegative");
		return value;
	}

	private static int residue(String s){
		require(s.length()==1, "Bad amino-acid symbol: "+s);
		final char c=s.charAt(0);
		final int code=c<AminoAcid.acidToNumber.length ? AminoAcid.acidToNumber[c] : -1;
		require(code>=0 && code<20 && AminoAcid.numberToAcid[code]==c, "Bad amino-acid symbol: "+s);
		return code;
	}

	private static int dimer(String s){
		require(s.length()==2, "Bad dimer symbol: "+s);
		return 20*residue(s.substring(0, 1))+residue(s.substring(1, 2));
	}

	private static boolean equals(byte[] line, byte[] expected){
		return Arrays.equals(line, expected);
	}

	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		writer.start();
		return writer;
	}

	private static void write(Path path, String text){
		try{Files.write(path, text.getBytes(StandardCharsets.US_ASCII), java.nio.file.StandardOpenOption.CREATE_NEW);}
		catch(java.io.IOException failure){throw new RuntimeException("Output failed: "+path, failure);}
	}

	private static String required(HashMap<String,String> options, String key){
		final String value=options.get(key);
		require(value!=null && value.length()>0, "Missing option: "+key);
		return value;
	}

	private static void require(boolean condition, String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static final class MonomerTable{
		MonomerTable(String[] bands_){
			bands=bands_;
			counts=new long[bands.length][20][20];
			alignments=new long[bands.length];
		}

		static MonomerTable add(MonomerTable result, MonomerTable source, long weight){
			if(result==null){result=new MonomerTable(source.bands);}
			require(Arrays.equals(result.bands, source.bands), "Mismatched monomer bands");
			for(int b=0; b<result.bands.length; b++){
				result.alignments[b]=Math.addExact(result.alignments[b], Math.multiplyExact(source.alignments[b], weight));
				for(int q=0; q<20; q++){for(int r=0; r<20; r++){
					result.counts[b][q][r]=Math.addExact(result.counts[b][q][r],
						Math.multiplyExact(source.counts[b][q][r], weight));
				}}
			}
			return result;
		}

		final String[] bands;
		final long[][][] counts;
		final long[] alignments;
	}

	private static final class DimerTable{
		DimerTable(String[] bands_){
			bands=bands_;
			ordered=new long[bands.length][DIMERS*DIMERS];
		}

		static DimerTable add(DimerTable result, DimerTable source, long weight){
			if(result==null){result=new DimerTable(source.bands);}
			require(Arrays.equals(result.bands, source.bands), "Mismatched dimer bands");
			for(int b=0; b<result.bands.length; b++){
				for(int i=0; i<DIMERS*DIMERS; i++){
					result.ordered[b][i]=Math.addExact(result.ordered[b][i], Math.multiplyExact(source.ordered[b][i], weight));
				}
			}
			return result;
		}

		final String[] bands;
		final long[][] ordered;
	}

	private static final int DIMERS=FamilyDimerSubstitutionTable.DIMERS;
	private static final String MONOMER_HEADER_TEXT=
		"band\tquery_aa\tconsensus_aa\tcount\tquery_monomer\tconsensus_monomer\talignments\taligned_pairs";
	private static final String DIMER_HEADER_TEXT="band\tquery_dimer\tconsensus_dimer\tcount";
	private static final byte[] MONOMER_HEADER=MONOMER_HEADER_TEXT.getBytes(StandardCharsets.US_ASCII);
	private static final byte[] DIMER_HEADER=DIMER_HEADER_TEXT.getBytes(StandardCharsets.US_ASCII);
	private static final HashMap<String,String> OPTIONS=new HashMap<String,String>();
	static{
		for(String key : new String[]{"monomers", "dimers", "weights", "out"}){
			OPTIONS.put(key, key);
		}
	}
}
