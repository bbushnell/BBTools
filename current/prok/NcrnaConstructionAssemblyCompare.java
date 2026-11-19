package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Positive whole-assembly identity evidence for explicitly nominated pairs.
 * Compares complete contig multisets by bytes, ignoring names, order, case and
 * each contig's orientation. Multiplicity and DNA IUPAC symbols are preserved.
 * A nonmatch remains unresolved: circular rotation, alternate fragmentation,
 * assembly versions and biological relatedness are outside this exact criterion.
 * Neither short-gene identity nor hashes establish an assembly alias here.
 * @author Ganyu
 */
public final class NcrnaConstructionAssemblyCompare{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssemblyCompare.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			int expected=-1;
			long maxBases=100000000;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("expectedpairs")){expected=Parse.parseIntKMG(b);}
				else if(a.equals("maxassemblybases")){maxBases=Parse.parseKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown assembly-comparison parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use -ea in=pairs.tsv out=fresh_prefix expectedpairs=N without sampling or overwrite");
			process(parser.in1, parser.out1, expected, maxBases);
			pp.outstream.println("NCRNA_ASSEMBLY_SEQUENCE_COMPARISON_COMPLETE_POSITIVE_IDENTITY_ONLY");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Retains every nominated pair; only positive exact matches establish aliases. */
	public static void process(final String input, final String output, final int expected, final long maxBases){
		require(input!=null && output!=null && expected>=0 && expected<=1000000 && maxBases>0,
			"Comparison needs a bounded declared pair count and positive per-assembly base limit");
		checkEnvironment();
		RrnaResourceIO.fresh(output, ".pairs.tsv.gz", ".summary.tsv");
		final ArrayList<Pair> pairs=loadPairs(input, expected, maxBases);
		Assembly left=null;
		Assembly right=null;
		long matched=0;
		long reads=0;
		long bases=0;
		final ByteBuilder buffer=new ByteBuilder();
		try(Output rows=new Output(output+".pairs.tsv.gz"); Output summary=new Output(output+".summary.tsv")){
			rows.writer.println(HEADER+"\tsequence_identity\tcriterion");
			for(final Pair pair:pairs){
				if(left==null || !left.path.equals(pair.left)){
					left=null;//Release the preceding payload before reading its replacement.
					left=read(pair.left, pair.bases, pair.contigs);
					reads++;
					bases+=pair.bases;
				}
				if(right==null || !right.path.equals(pair.right)){
					right=null;
					right=read(pair.right, pair.bases, pair.contigs);
					reads++;
					bases+=pair.bases;
				}
				final boolean same=equal(left.contigs, right.contigs);
				if(same){matched++;}
				rows.writer.print(buffer.clear().append(pair.row).tab().append(same ? MATCH : UNRESOLVED)
					.tab().append(CRITERION).nl());
			}
			summary.writer.println("pairs\texact_matches\tunresolved_nonmatches\tassembly_reads\tbases_read\tstatus");
			summary.writer.print(buffer.clear().append(pairs.size()).tab().append(matched).tab()
				.append(pairs.size()-matched).tab().append(reads).tab().append(bases)
				.append("\tPOSITIVE_IDENTITY_ONLY_NOT_FROZEN_SPLIT\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Pair Provenance       ----------------*/
	/*--------------------------------------------------------------*/

	private static ArrayList<Pair> loadPairs(final String path, final int expected, final long maxBases){
		final ArrayList<Pair> pairs=new ArrayList<Pair>();
		final ObjectMap<String, Pair> identities=new ObjectMap<String, Pair>(String.class, Pair.class);
		final ObjectMap<String, Pair> paths=new ObjectMap<String, Pair>(String.class, Pair.class);
		final ObjectMap<String, String> inputs=new ObjectMap<String, String>(String.class, String.class);
		final LineParser1 parser=new LineParser1('\t');
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==7 && parser.termEquals(BASIS, 6), "Pairs must come from the dimension-candidate audit: "+path);
				final String eid=utf8(parser, row, 0);
				final String did=utf8(parser, row, 1);
				final Pair pair=new Pair(row, utf8(parser, row, 2), utf8(parser, row, 3), parser.parseLong(4), parser.parseInt(5));
				require(!eid.isEmpty() && !did.isEmpty() && pair.bases>0 && pair.bases<=maxBases
					&& pair.contigs>0 && pair.contigs<=pair.bases,
					"Each pair needs two cohort-qualified IDs and bounded positive complete-assembly dimensions");
				final String key=eid+'\t'+did;
				require(identities.get(key)==null, "Duplicate evaluation/discovery pair: "+key);
				identities.put(key, pair);
				bind(inputs, "evaluation\t"+eid, pair.left);
				bind(inputs, "discovery\t"+did, pair.right);
				bindPath(paths, pair.left, pair);
				bindPath(paths, pair.right, pair);
				pairs.add(pair);
				require(pairs.size()<=expected, "Pair input exceeds its declared complete count");
			}
		}
		require(pairs.size()==expected, "Pair input is missing declared comparisons");
		return pairs;
	}

	private static void bind(final ObjectMap<String, String> inputs, final String key, final String path){
		final String old=inputs.get(key);
		require(old==null || old.equals(path), "One cohort-qualified input cannot change its FASTA path: "+key);
		inputs.put(key, path);
	}

	private static void bindPath(final ObjectMap<String, Pair> paths, final String path, final Pair pair){
		final File file=new File(path);
		require(file.isAbsolute() && file.isFile(), "A nominated assembly needs an existing absolute FASTA path: "+path);
		final Pair old=paths.get(path);
		require(old==null || old.bases==pair.bases && old.contigs==pair.contigs,
			"One assembly path cannot change dimensions between nominated pairs: "+path);
		paths.put(path, pair);
	}

	private static String utf8(final LineParser1 parser, final byte[] row, final int column){
		parser.setBounds(column);
		final String value=new String(row, parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 must not silently change input identity");
		return value;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Exact Sequences        ----------------*/
	/*--------------------------------------------------------------*/

	private static void checkEnvironment(){
		require(!Shared.AMINO_IN && !Read.IUPAC_TO_N && !Read.U_TO_T && !Read.DOT_DASH_X_TO_N
			&& !Read.LOWER_CASE_TO_N && Read.JUNK_MODE==Read.CRASH_JUNK,
			"Whole-assembly identity forbids amino mode, masking or lossy nucleotide conversion");
	}

	private static Assembly read(final String path, final long expectedBases, final int expectedContigs){
		final File file=new File(path);
		final long bytes=file.length();
		final long modified=file.lastModified();
		final FileFormat format=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		require(format.fasta() && !format.interleaved() && !format.amino(), "Assembly input must be unpaired nucleotide FASTA: "+path);
		final ArrayList<byte[]> contigs=new ArrayList<byte[]>();
		long count=0;
		//Zero reader threads selects a synchronous native FASTA implementation;
		//a validation exception cannot strand an asynchronous parsing producer.
		try(SequenceInput in=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 0))){
			in.stream.start();
			for(ListNum<Read> batch=in.stream.nextList(); batch!=null; batch=in.stream.nextList()){
				for(final Read read:batch.list){
					require(read.mate==null && read.bases!=null && read.length()>0, "Empty or paired assembly record: "+path);
					count+=read.length();
					require(count<=expectedBases && contigs.size()<expectedContigs, "Assembly exceeds its audited dimensions: "+path);
					canonical(read.bases);
					contigs.add(read.bases);
				}
				in.stream.returnList(batch);
			}
		}
		require(count==expectedBases && contigs.size()==expectedContigs, "Complete assembly differs from its audited dimensions: "+path);
		require(file.length()==bytes && file.lastModified()==modified, "Assembly file changed during comparison: "+path);
		Collections.sort(contigs, ORDER);
		return new Assembly(path, contigs);
	}

	/** Owns these native Read bases. All accepted symbols have involutive DNA complements. */
	private static void canonical(final byte[] sequence){
		assert(sequence!=null && sequence.length>0) : "Assembly records are checked before orientation normalization";
		for(int i=0; i<sequence.length; i++){
			final int base=sequence[i]&255;
			final byte upper=(byte)(base>='a' && base<='z' ? base-32 : base);
			require(base<128 && "ACGTRYSWKMBDHVN".indexOf(upper)>=0, "Only complete DNA IUPAC symbols are supported; byte="+base);
			sequence[i]=upper;
		}
		for(int i=0; i<sequence.length; i++){
			final int difference=sequence[i]-AminoAcid.baseToComplementExtended[sequence[sequence.length-1-i]];
			if(difference<0){return;}
			if(difference>0){
				AminoAcid.reverseComplementBasesInPlace(sequence);
				return;
			}
		}
	}

	private static boolean equal(final ArrayList<byte[]> a, final ArrayList<byte[]> b){
		if(a.size()!=b.size()){return false;}
		for(int i=0; i<a.size(); i++){
			if(!Arrays.equals(a.get(i), b.get(i))){return false;}
		}
		return true;
	}

	private static final Comparator<byte[]> ORDER=new Comparator<byte[]>(){
		@Override
		public int compare(final byte[] a, final byte[] b){
			if(a.length!=b.length){return Integer.compare(a.length, b.length);}
			for(int i=0; i<a.length; i++){
				if(a[i]!=b[i]){return a[i]-b[i];}
			}
			return 0;
		}
	};

	private static final class Pair{
		Pair(final byte[] row_, final String left_, final String right_, final long bases_, final int contigs_){
			row=row_;
			left=left_;
			right=right_;
			bases=bases_;
			contigs=contigs_;
		}
		final byte[] row;
		final String left, right;
		final long bases;
		final int contigs;
	}
	private static final class Assembly{
		Assembly(final String path_, final ArrayList<byte[]> contigs_){
			path=path_;
			contigs=contigs_;
		}
		final String path;
		final ArrayList<byte[]> contigs;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	public static final String HEADER="evaluation_input\tdiscovery_input\tevaluation_fasta\tdiscovery_fasta\tbases\tcontigs\tbasis";
	public static final String BASIS="DIMENSION_COMPATIBILITY_NOT_SEQUENCE_IDENTITY";
	public static final String MATCH="EXACT_ASSEMBLY_SEQUENCE_ALIAS";
	public static final String UNRESOLVED="NONMATCH_BIOLOGICAL_IDENTITY_UNRESOLVED";
	public static final String CRITERION="DNA_IUPAC_CONTIG_MULTISET_CASE_ORDER_ORIENTATION_NORMALIZED";
}
