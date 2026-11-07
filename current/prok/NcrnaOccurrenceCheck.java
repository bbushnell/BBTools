package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.FileFormat;
import map.IntHashMap;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaObservationData.Symbols;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/** Checks actual extracted contexts/cores against their complete occurrence ledger.
 * Native sequence readers finish before validation, so a metadata mismatch cannot
 * strand a producer behind the checker. Both short-RNA FASTAs are retained; memory
 * is proportional to their total bases. No genome reread or biological filter.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOccurrenceCheck{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaOccurrenceCheck tool=new NcrnaOccurrenceCheck(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	public NcrnaOccurrenceCheck(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			final Parser parser=new Parser();
			for(final String arg:args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}
				if(a.equals("contexts")){contexts=b;}
				else if(a.equals("cores")){cores=b;}
				else if(a.equals("expected")){expected=Parse.parseIntKMG(b);}
				else if(parser.parse(arg, a, b)){
					// Native shared options.
				}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			in=parser.in1;
			out=parser.out1;
			require(in!=null && contexts!=null && cores!=null && out!=null && expected>0,
				"Required: in=occurrence_loci.tsv contexts=contexts.fa cores=cores.fa out=fresh_summary.tsv expected=N");
			require(!parser.append && parser.maxReads<0 && parser.samplerate==1 && parser.in2==null && parser.out2==null,
				"Check the complete unpaired population without sampling or append");
			require(Tools.testOutputFiles(false, false, false, out), "A fresh writable summary is required");
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	private void process(final Timer timer){
		assert(rows.isEmpty()) : "A checker owns one complete declared occurrence population";
		loadLedger();
		final boolean trim=Shared.TRIM_READ_DESCRIPTION, upper=Read.TO_UPPER_CASE;
		final boolean validate=Read.VALIDATE_IN_CONSTRUCTOR;
		Shared.TRIM_READ_DESCRIPTION=false;
		Read.TO_UPPER_CASE=false;
		// Read constructor validation can repair symbols in place. Parity compares
		// actual bytes instead; native I/O error checks remain active in toReads.
		Read.VALIDATE_IN_CONSTRUCTOR=false;
		try{
			checkContexts(ReadInputStream.toReads(contexts, FileFormat.FASTA, -1));
			checkCores(ReadInputStream.toReads(cores, FileFormat.FASTA, -1));
		}finally{
			Shared.TRIM_READ_DESCRIPTION=trim;
			Read.TO_UPPER_CASE=upper;
			Read.VALIDATE_IN_CONSTRUCTOR=validate;
		}
		try(Output output=new Output(out)){
			output.writer.println("occurrences\tcontexts\tcores\tcontext_bases\tcore_bases\tcontigs_with_gc");
			output.writer.print(buffer.clear().append(rows.size()).tab().append(rows.size()).tab().append(rows.size())
				.tab().append(contextBases).tab().append(coreBases).tab().append(gcByContig.size()).nl());
		}
		timer.stop();
		outstream.println("NCRNA_OCCURRENCE_SEQUENCE_PARITY_PASS occurrences="+rows.size()+" time="+timer);
	}

	/** Retains one ID and geometry row per occurrence, with byte-interned sources. */
	private void loadLedger(){
		assert(rows.isEmpty()) : "Ledger IDs are established before sequence joins";
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, NcrnaLearningPopulation.OCCURRENCE_HEADER, in);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==15, "Expected complete occurrence metadata");
				final String id=parser.parseString(0);
				require(!id.isEmpty() && !byId.contains(id), "Record IDs must be nonempty and unique");
				parser.setBounds(1);
				require(parser.b()>parser.a() && id.length()==parser.b()-parser.a()+3 && id.startsWith("ID="), "Record and occurrence IDs disagree");
				for(int i=parser.a(); i<parser.b(); i++){
					require(id.charAt(i-parser.a()+3)==(line[i]&255), "Record and occurrence IDs disagree");
				}
				parser.setBounds(2);
				final int source=symbols.intern(line, parser.a(), parser.b());
				buffer.clear();
				for(final int field:SOURCE_FIELDS){
					parser.setBounds(field);
					require(parser.b()>parser.a(), "Source, assembly and contig identities must not be empty");
					buffer.append(line, parser.a(), parser.b()-parser.a()).tab();
				}
				final int contig=symbols.intern(buffer.array, 0, buffer.length());
				final int start=number(9), stop=number(10), length=number(11), contigLength=number(12);
				final int left=number(13), right=number(14);
				final boolean minus=parser.termEquals('-', 8);
				require(minus || parser.termEquals('+', 8), "Explicit occurrence strand required");
				require(start>0 && stop>=start && stop<=contigLength && (long)stop-start+1==length,
					"Ledger coordinates and length disagree");
				require(left<=(minus ? contigLength-stop : start-1) && right<=(minus ? start-1 : contigLength-stop),
					"Sense flanks exceed the declared contig");
				byId.put(id, rows.size());
				rows.add(new Row(id, source, contig, length, left, right));
			}
		}
		require(rows.size()==expected, "Ledger does not contain the complete declared population");
	}

	private void checkContexts(final ArrayList<Read> reads){
		require(reads.size()==expected, "Context FASTA count differs from the complete ledger");
		for(final Read read:reads){
			final Row row=lookup(read);
			if(row.context!=null){throw new IllegalArgumentException("Duplicated context record: "+row.id);}
			if(read.bases==null || (long)read.length()!=(long)row.length+row.left+row.right){
				throw new IllegalArgumentException("Actual context length disagrees with declared core/flanks: "+row.id);
			}
			checkNamedHeader(read.id, row, true);
			row.context=read;
			contextBases+=read.length();
		}
		for(final Row row:rows){require(row.context!=null, "Missing context record");}
	}

	private void checkCores(final ArrayList<Read> reads){
		require(reads.size()==expected, "Core FASTA count differs from the complete ledger");
		for(final Read read:reads){
			final Row row=lookup(read);
			if(row.coreSeen || read.bases==null || read.length()!=row.length){
				throw new IllegalArgumentException("Duplicate core or incorrect exact-core length: "+row.id);
			}
			checkNamedHeader(read.id, row, false);
			assert(row.context!=null) : "The complete context pass precedes any core comparison";
			for(int i=0; i<read.length(); i++){
				if(read.bases[i]!=row.context.bases[row.left+i]){
					throw new IllegalArgumentException("Core differs from its actual context at "+row.id+" base "+i);
				}
			}
			row.coreSeen=true;
			coreBases+=read.length();
		}
		for(final Row row:rows){require(row.coreSeen, "Missing exact core");}
	}

	/** Actual header values must agree with metadata; GC also agrees across a source contig.
	 * GC is the caller-compatible value supplied by CutGff, not an independent genome recount.
	 */
	private void checkHeader(final String text, final Row row, final boolean context){
		assert(text.startsWith(row.id)) : "lookup must bind the exact record ID before header scanning starts at its length";
		int seen=0;
		for(int from=row.id.length(); from<text.length();){
			while(from<text.length() && Character.isWhitespace(text.charAt(from))){from++;}
			int to=from;
			while(to<text.length() && !Character.isWhitespace(text.charAt(to))){to++;}
			int flag=0;
			if(text.startsWith("source=", from)){
				flag=1;
				final byte[] expected=symbols.bytes(row.source);
				require(to-from-7==expected.length, "Context/core source differs from its ledger");
				for(int i=0; i<expected.length; i++){require(text.charAt(from+7+i)==(expected[i]&255), "Source identity mismatch");}
			}else if(text.startsWith("lflank=", from)){
				flag=2;
				require(NcrnaNumbers.unsignedInt(text, from+7, to)==row.left, "Actual left flank disagrees with metadata");
			}else if(text.startsWith("rflank=", from)){
				flag=4;
				require(NcrnaNumbers.unsignedInt(text, from+7, to)==row.right, "Actual right flank disagrees with metadata");
			}else if(text.startsWith("contig_gc=", from)){
				flag=8;
				final float gc=Float.parseFloat(text.substring(from+10, to));
				require(Float.isFinite(gc) && gc>=0 && gc<=1, "Actual contig GC must be finite and in [0,1]");
				final int bits=gc==0 ? 0 : Float.floatToIntBits(gc), prior=gcByContig.get(row.contig);
				require(prior<0 || prior==bits, "One source contig has inconsistent caller GC values");
				gcByContig.put(row.contig, bits);
			}
			require((seen&flag)==0, "Duplicated source/flank/GC header field");
			seen|=flag;
			from=to;
		}
		require(context ? seen==15 : seen==1, "Context needs source/flanks/GC; exact core needs source only");
	}

	/** Adds the exact failing identity without allocating diagnostic text for valid records. */
	private void checkNamedHeader(final String text, final Row row, final boolean context){
		try{checkHeader(text, row, context);}
		catch(IllegalArgumentException failure){
			throw new IllegalArgumentException(row.id+": "+failure.getMessage(), failure);
		}
	}

	private Row lookup(final Read read){
		require(read!=null && read.mate==null, "Expected unpaired occurrence sequences");
		final String id=NcrnaSupport.recordName(read.id);
		final int index=byId.get(id);
		if(index<0){throw new IllegalArgumentException("Sequence ID is absent from the complete occurrence ledger: "+id);}
		return rows.get(index);
	}
	private int number(final int field){
		parser.setBounds(field);
		return NcrnaNumbers.unsignedInt(parser.line(), parser.a(), parser.b());
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static final class Row{
		Row(String id_, int source_, int contig_, int length_, int left_, int right_){
			id=id_;
			source=source_;
			contig=contig_;
			length=length_;
			left=left_;
			right=right_;
		}
		final String id;
		final int source, contig, length, left, right;
		Read context;
		boolean coreSeen;
	}
	private String in, contexts, cores, out;
	private int expected=-1;
	private long contextBases, coreBases;
	private PrintStream outstream;
	private final ArrayList<Row> rows=new ArrayList<Row>();
	private final ObjectIntMap<String> byId=new ObjectIntMap<String>(String.class);
	private final Symbols symbols=new Symbols();
	private final IntHashMap gcByContig=new IntHashMap();
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder buffer=new ByteBuilder();
	private static final int[] SOURCE_FIELDS={2, 3, 7};
}
