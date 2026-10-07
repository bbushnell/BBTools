package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Inventories declared CDS translation tables in one flat paired GFF/FASTA directory.
 * Missing declarations remain table 0; no biological code is inferred from taxonomy.
 * Counts are GFF CDS rows, not deduplicated protein identifiers. Input is never modified.
 * @author Keqing
 */
public final class TranslationTableCensus {

	public static void main(String[] args){
		final TranslationTableCensus census=new TranslationTableCensus(args);
		try{census.process();}finally{Shared.closeStream(census.outstream);}
	}

	/** Parses standard BBTools options plus a flat input directory. Output is written after all reads succeed. */
	private TranslationTableCensus(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg : pp.args){
			final String[] pair=arg.split("=", 2);
			final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length==2 ? pair[1] : null;
			if(key.equals("indir")){directory=value;}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		if(directory==null || parser.out1==null){throw new IllegalArgumentException("Require indir=DIR and out=TABLE.tsv");}
		if(parser.append){throw new IllegalArgumentException("A census is a complete inventory; append is unsupported");}
		outputFile=new File(parser.out1).getAbsoluteFile();
		if(!Tools.testOutputFiles(parser.overwrite, false, false, parser.out1)){
			throw new IllegalArgumentException("Cannot write census: "+parser.out1);
		}
		final File dir=new File(directory).getAbsoluteFile();
		files=dir.listFiles((parent, name)->name.endsWith(".gff") || name.endsWith(".gff.gz"));
		if(files==null || files.length==0){throw new IllegalArgumentException("No GFF files in readable directory: "+dir);}
		Arrays.sort(files);
		results=new Counts[files.length];
		pairs=new File[files.length];
		for(File file : files){
			checkPath(file);
			if(!file.isFile() || !file.canRead()){throw new IllegalArgumentException("Unreadable GFF: "+file);}
			checkOutputAlias(file);
		}
		ffout=FileFormat.testOutput(parser.out1, FileFormat.TXT, null, false, parser.overwrite, false, false);
	}

	/** Independent file workers retain only small counters. No successful output can hide a failed worker. */
	private void process(){
		final Worker[] workers=new Worker[Math.min(files.length, Math.max(1, Shared.threads()))];
		for(int i=0; i<workers.length; i++){workers[i]=new Worker(); workers[i].start();}
		for(Worker worker : workers){
			try{worker.join();}catch(InterruptedException e){Thread.currentThread().interrupt(); throw new RuntimeException("Census interrupted", e);}
		}
		for(Worker worker : workers){
			if(worker.failure!=null){throw new RuntimeException("Census failed; no inventory was written", worker.failure);}
		}
		final ByteStreamWriter writer=ByteStreamWriter.makeBSW(ffout);
		final ByteBuilder row=new ByteBuilder();
		long total=0;
		try{
			writer.print("gff\tfasta\tpair_status\tcode_status\ttable\tcds_rows\tcds_total\tmissing_rows\texception_rows\n");
			for(int i=0; i<files.length; i++){
				final Counts count=results[i];
				assert(count!=null) : "Joined workers must populate every inventory slot; missing "+files[i];
				total+=count.total;
				final File fasta=pairs[i];
				final String state=count.total==0 ? "NO_CDS" : count.missing>0 ? "MISSING" : count.codes.size>1 ? "MIXED" : "UNIFORM";
				if(count.codes.size==0){count.codes.add(0); count.rows.add(0);}
				long sum=0;
				for(int j=0; j<count.codes.size; j++){
					sum+=count.rows.get(j);
					row.clear().append(files[i].toString()).tab().append(fasta==null ? "." : fasta.toString()).tab()
						.append(fasta==null ? "MISSING" : "OK").tab().append(state).tab().append(count.codes.get(j)).tab()
						.append(count.rows.get(j)).tab().append(count.total).tab().append(count.missing).tab().append(count.exceptions).nl();
					writer.print(row);
				}
				assert(sum==count.total) : "Every CDS must belong to one declared/missing code bucket: "+files[i];
			}
		}finally{
			if(writer.poisonAndWait()){throw new RuntimeException("Census output failed: "+ffout.name());}
		}
		outstream.println("CENSUS_COMPLETE files="+files.length+" cds_rows="+total);
	}

	/** Input identity includes symlinks and hard links; ow=t must never authorize overwriting a corpus file. */
	private void checkOutputAlias(final File input){
		try{
			if(outputFile.exists() && Files.isSameFile(input.toPath(), outputFile.toPath())){
				throw new IllegalArgumentException("Census output aliases an input: "+input);
			}
		}catch(IOException e){throw new IllegalArgumentException("Cannot check input/output identity: "+input, e);}
	}

	/** Checks numbered-code training truth before a trainer reopens the GFF.
	 * A regular readable file is required because validation and training are separate reads. */
	static void requireDeclaredCode(final String path, final int table){
		if(table<1){throw new IllegalArgumentException("A numbered translation table must be positive: "+table);}
		final File file=new File(path);
		if(!file.isFile() || !file.canRead()){
			throw new IllegalArgumentException("Numbered-table training requires a readable, rereadable GFF file: "+path);
		}
		final Counts counts=scan(file);
		if(counts.total==0 || counts.missing!=0 || counts.codes.size!=1 || counts.codes.get(0)!=table){
			throw new IllegalArgumentException("Training requires uniform declared transl_table="+table
				+" on every CDS: "+path+"; CDS="+counts.total+", missing="+counts.missing+", distinct="+counts.codes.size);
		}
	}

	/** Reads one GFF through ByteFile; embedded FASTA terminates annotation parsing. */
	static Counts scan(final File file){
		final ByteFile reader=ByteFile.makeByteFile1(file.toString(), false);
		final LineParser1 fields=new LineParser1('\t');
		final Counts counts=new Counts();
		long number=0;
		try{
			for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
				number++;
				if(line.length==0){continue;}
				if(line[0]=='#'){
					if(matches(line, 0, line.length, "##FASTA")){break;}
					continue;
				}
				fields.set(line);
				if(fields.terms()!=9){throw new IllegalArgumentException("Expected 9 GFF columns: "+file+":"+number);}
				if(!fields.termEquals("CDS", 2)){continue;}
				fields.setBounds(8);
				try{counts.add(line, fields.a(), fields.b());}
				catch(IllegalArgumentException e){throw new IllegalArgumentException(file+":"+number+": "+e.getMessage(), e);}
			}
		}finally{
			if(reader.close()){throw new RuntimeException("GFF read failed: "+file+":"+number);}
		}
		return counts;
	}

	/** Resolves only exact basename pairs. Multiple candidate FASTAs are ambiguous and fail loudly. */
	static File pairedFasta(final File gff){
		final String name=gff.toString();
		final String prefix=name.substring(0, name.length()-(name.endsWith(".gz") ? 7 : 4));
		File found=null;
		for(String suffix : new String[]{".fna.gz", ".fna", ".fa.gz", ".fa", ".fasta.gz", ".fasta"}){
			final File candidate=new File(prefix+suffix);
			if(candidate.exists()){
				if(found!=null){throw new IllegalArgumentException("Ambiguous FASTA pair: "+found+" and "+candidate);}
				if(!candidate.isFile() || !candidate.canRead()){throw new IllegalArgumentException("Unreadable FASTA pair: "+candidate);}
				checkPath(candidate);
				found=candidate;
			}
		}
		return found;
	}

	/** Rejects paths that cannot be represented unambiguously in the output TSV. */
	private static void checkPath(final File file){
		final String path=file.toString();
		if(path.indexOf('\t')>=0 || path.indexOf('\n')>=0 || path.indexOf('\r')>=0){
			throw new IllegalArgumentException("Control character in inventory path: "+file);
		}
	}

	/** Exact ASCII attribute-key comparison, avoiding substring and per-record String parsing. */
	private static boolean matches(final byte[] line, final int start, final int end, final String value){
		if(end-start!=value.length()){return false;}
		for(int i=0; i<value.length(); i++){if(line[start+i]!=value.charAt(i)){return false;}}
		return true;
	}

	/** Small sparse primitive histogram. Zero denotes missing metadata, never code 11. */
	static final class Counts {
		void add(final byte[] line, final int start, final int end){
			int code=0;
			boolean seen=false, exception=false;
			for(int a=start; a<end;){
				int b=a;
				while(b<end && line[b]!=';'){b++;}
				int eq=a;
				while(eq<b && line[eq]!='='){eq++;}
				if(matches(line, a, eq, "transl_table")){
					if(seen){throw new IllegalArgumentException("Duplicate transl_table attribute");}
					seen=true;
					if(eq+1>=b){throw new IllegalArgumentException("Empty transl_table attribute");}
					for(int p=eq+1; p<b; p++){
						final int digit=line[p]-'0';
						if(digit<0 || digit>9 || code>(Integer.MAX_VALUE-digit)/10){
							throw new IllegalArgumentException("Invalid positive integer transl_table attribute");
						}
						code=10*code+digit;
					}
					if(code<1){throw new IllegalArgumentException("transl_table must be positive");}
				}else if(matches(line, a, eq, "transl_except") || matches(line, a, eq, "exception")){exception=true;}
				a=b+1;
			}
			total++;
			if(!seen){missing++;}
			if(exception){exceptions++;}
			int index=0;
			while(index<codes.size && codes.get(index)!=code){index++;}
			if(index==codes.size){codes.add(code); rows.add(0);}
			rows.increment(index);
			assert(total>=missing && total>=exceptions) : "Subset CDS counters exceed all CDS rows";
		}
		final IntList codes=new IntList();
		final LongList rows=new LongList();
		long total, missing, exceptions;
	}

	/** Each worker owns its parser and counters; slots are disjoint and consumed only after joins. */
	private final class Worker extends Thread {
		@Override
		public void run(){
			try{
				for(int i=next.getAndIncrement(); i<files.length; i=next.getAndIncrement()){
					pairs[i]=pairedFasta(files[i]);
					if(pairs[i]!=null){checkOutputAlias(pairs[i]);}
					results[i]=scan(files[i]);
					final int n=completed.incrementAndGet();
					if(n%100==0){outstream.println("CENSUS_PROGRESS files="+n+"/"+files.length);}
				}
			}catch(Throwable e){failure=e; next.set(files.length);}
		}
		Throwable failure;
	}

	private String directory;
	private final PrintStream outstream;
	private final File[] files;
	private final Counts[] results;
	private final File[] pairs;
	private final FileFormat ffout;
	private final File outputFile;
	private final AtomicInteger next=new AtomicInteger(), completed=new AtomicInteger();
}
