package test.jgi;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

import dna.Data;
import fileIO.ReadWrite;
import jgi.ReformatStreamer;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import processor.ReformatProcessor;
import shared.Shared;
import stream.FASTQ;
import stream.Read;
import stream.ReadStreamWriter;
import stream.SamLine;
import stream.SamWriter;
import stream.Writer;
import structures.ListNum;
import tracker.ReadStats;

/** Baseline and strand-repair SAM splitting observations with literal strand/quality oracles.
 * Covers the real SAM-to-SAM route, direct production Read-batch processing plus
 * real SAM serialization, and SAM-to-FASTQ conversion. Baseline records inherited
 * reverse flags; fixed requires load-flipped fragments to inherit forward strand.
 * Baseline retains the historical unsplit SAM route; fixed also requires the
 * subsequently accepted SAM-to-SAM splitting. Short-record metadata and downstream rcomp remain.
 * Invoke through the dev launcher in a fresh JVM with a new output directory.
 * @author Sucrose
 */
public final class ReformatSamBreakBaselineTest{

	/** Requires phase=baseline or fixed and outdir; generates tiny unpaired primary records. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if("phase".equals(pair[0])){phase=pair[1];}
			else if("outdir".equals(pair[0])){directory=pair[1];}
			else{throw new IllegalArgumentException(arg);}
		}
		check("baseline".equals(phase) || "fixed".equals(phase), "Set phase=baseline or phase=fixed.");
		final boolean fixed="fixed".equals(phase);
		check(directory!=null, "Set outdir to a new directory with an existing parent.");
		check(Data.GENOME_BUILD<0 && Data.scaffoldLocs==null, "Baseline requires a fresh JVM without a loaded reference.");
		check(!ReadStats.collectingStats(), "Baseline requires default disabled statistics in its fresh JVM.");
		final Path launcher=Paths.get("reformat.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "Missing actual BBTools launcher: "+launcher);
		final Path root=Files.createDirectory(Paths.get(directory).toAbsolutePath());
		final Case[] cases=cases();
		final Path input=root.resolve("input.sam");
		final ArrayList<String> records=new ArrayList<String>();
		for(final Case item : cases){records.add(item.inputLine());}
		final ArrayList<String> sam=new ArrayList<String>(Arrays.asList("@HD\tVN:1.6\tSO:unknown", "@SQ\tSN:chr1\tLN:1000"));
		sam.addAll(records);
		Files.write(input, sam, StandardCharsets.US_ASCII);
		Files.write(root.resolve("expected-sam-route-records.txt"), records, StandardCharsets.US_ASCII);
		for(final boolean flip : new boolean[]{true, false}){
			final Path work=Files.createDirectory(root.resolve("flip-"+flip));
			directReads(input, work, cases, flip, false, fixed);
			samRoute(launcher, input, work, records, cases, flip, fixed);
			fastqRoute(launcher, input, work, cases, flip);
		}
		//Retain direct/FASTQ controls above; the fixed SAM route also checks the later splitting policy.
		final Case x6Long=new Case("x6_mr_long", true, true, false, Integer.valueOf(1));
		final Case x6Short=new Case("x6_mr_short", true, true, true, Integer.valueOf(1));
		for(final boolean flip : new boolean[]{true, false}){
			final Path rcWork=Files.createDirectory(root.resolve("rcomp-flip-"+flip));
			directReads(writeCaseInput(rcWork, cases), rcWork, cases, flip, true, fixed);
			final Case[] x6Cases={x6Long, x6Short};
			final Path x6Work=Files.createDirectory(root.resolve("x6-flip-"+flip));
			directReads(writeCaseInput(x6Work, x6Cases), x6Work, x6Cases, flip, false, fixed);
			final Case[] x6RcCases={x6Long, x6Short};
			final Path x6RcWork=Files.createDirectory(root.resolve("x6-rcomp-flip-"+flip));
			directReads(writeCaseInput(x6RcWork, x6RcCases), x6RcWork, x6RcCases, flip, true, fixed);
		}
		System.out.println("REFORMAT_SAM_BREAK_"+(fixed ? "FIXED" : "BASELINE")+"_PASS checks="+checks);
	}

	/** Retains the exact extra input records supplied to the direct production Read route. */
	private static Path writeCaseInput(final Path work, final Case[] cases) throws Exception{
		final ArrayList<String> lines=new ArrayList<String>(Arrays.asList("@HD\tVN:1.6\tSO:unknown", "@SQ\tSN:chr1\tLN:1000"));
		for(final Case item : cases){lines.add(item.inputLine());}
		final Path input=work.resolve("extension-input.sam");
		Files.write(input, lines, StandardCharsets.US_ASCII);
		return input;
	}
	/** Extra rcomp cases retain current short-record behavior, including its stale attached SAM flags. */
	private static void directReads(final Path input, final Path work, final Case[] cases,
			final boolean flip, final boolean rcomp, final boolean fixed) throws Exception{
		final SavedStatics saved=new SavedStatics(SamLine.class, Read.class, ReadStreamWriter.class,
				ReadWrite.class, ReadStats.class, FASTQ.class, Parser.class, Shared.class);
		try{
			configureSam(flip);
			final ReformatStreamer driver=new ReformatStreamer(new String[]{"in="+input,
					"out="+work.resolve("constructor-only-unused.sam"), "wt=1", "tin=1", "tout=1",
					"int=f", "ow=f", "simd=f", "minlength=1", "breaklength=3", "flipsam="+flip, "rcomp="+rcomp});
			final Field field=ReformatStreamer.class.getDeclaredField("processor");
			field.setAccessible(true);
			final ReformatProcessor processor=(ReformatProcessor)field.get(driver);
			final ArrayList<Read> reads=new ArrayList<Read>();
			final ArrayList<Expected> expected=new ArrayList<Expected>();
			final LineParser1 parser=new LineParser1('\t');
			for(int i=0; i<cases.length; i++){
				final Case item=cases[i];
				final SamLine line=new SamLine(parser.set(item.inputLine().getBytes(StandardCharsets.US_ASCII)));
				final Read read=line.toRead(false);
				read.numericID=i;
				read.samline=line;
				line.obj=read;
				reads.add(read);
				final boolean reverse=item.reverseStored(flip);
				final String sequence=item.internalSequence(flip), quality=item.internalQuality(flip);
				System.out.println("OBS parsed flip="+flip+" rcomp="+rcomp+"\t"+details(read));
				check(sequence.equals(bases(read)) && quality.equals(qualities(read))
						&& read.mapped()==item.readMapped() && read.strand()==(item.readReverse() ? 1 : 0),
						"Unexpected parsed orientation for "+item.name+", flip="+flip+": "+details(read));
				if(item.shortRead){
					expected.add(new Expected(item.shortOutput(rcomp), read, line, i, false,
							item.readMapped(), item.readReverse()!=rcomp));
				}else{
					final String[] sequences=rcomp ? (reverse ? RC_REVERSE_PARTS : RC_FORWARD_PARTS)
							: (reverse ? REVERSE_PARTS : FORWARD_PARTS);
					final String[] qualities=rcomp ? (reverse ? RC_REVERSE_QUALS : RC_FORWARD_QUALS)
							: (reverse ? REVERSE_QUALS : FORWARD_QUALS);
					//This expectation uses literal input provenance, independently of production Read flags.
					final boolean fragmentReverse=(fixed && flip && item.mapped && item.reverse) ? false : item.readReverse();
					final int flag=(fragmentReverse!=rcomp) ? 20 : 4;
					for(int part=0; part<3; part++){
						expected.add(new Expected(unmappedLine(item.name+"_"+(part+1), flag,
								sequences[part], qualities[part]), read, line, i, true,
								false, fragmentReverse!=rcomp));
					}
				}
			}
			final RecordingWriter writer=new RecordingWriter();
			final Method method=ReformatStreamer.class.getDeclaredMethod("processReadList", ListNum.class,
					ReformatProcessor.class, Writer.class, Writer.class);
			method.setAccessible(true);
			try{method.invoke(driver, new ListNum<Read>(reads, 0), processor, writer, null);}
			catch(InvocationTargetException e){throw new AssertionError("Production Read-batch processing failed", e.getCause());}
			final ArrayList<SamLine> lines=SamWriter.toSamLines(writer.reads, true, false);
			final ArrayList<String> actualRows=new ArrayList<String>(), expectedRows=new ArrayList<String>();
			final ArrayList<String> observations=new ArrayList<String>();
			for(final Read read : writer.reads){observations.add(details(read));}
			for(final SamLine line : lines){actualRows.add(line.toString());}
			for(final Expected item : expected){expectedRows.add(item.row);}
			Files.write(work.resolve("direct.sam"), actualRows, StandardCharsets.US_ASCII);
			Files.write(work.resolve("expected-direct.sam"), expectedRows, StandardCharsets.US_ASCII);
			Files.write(work.resolve("direct-reads.tsv"), observations, StandardCharsets.US_ASCII);
			check(writer.ids.equals(Collections.singletonList(Long.valueOf(0)))
					&& writer.reads.size()==expected.size() && lines.size()==expected.size(),
					"Expected one batch containing"+expected.size()+" records; got "+writer.reads.size()+"/"+lines.size());
			for(int i=0; i<expected.size(); i++){
				final Expected item=expected.get(i);
				final Read read=writer.reads.get(i);
				check(read.mapped()==item.mapped && read.strand()==(item.reverse ? 1 : 0)
						&& (item.split ? read!=item.original && read.samline==null
						&& read.chrom==-1 && read.start==-1 && read.stop==-1
						: read==item.original && read.samline==item.line),
						"Read/attachment state differed at index "+i+": "+details(read));
				check(read.numericID==item.numericID, "Numeric ID changed at index "+i+": "+read.numericID);
				check(item.row.equals(actualRows.get(i)), "Serialized row differs at index "+i
						+"; expected="+item.row+"; actual="+actualRows.get(i));
			}
		}finally{saved.restore();}
	}

	/** Baseline retains eight unsplit rows; fixed checks complete split rows from literal input provenance. */
	private static void samRoute(final Path launcher, final Path input, final Path work,
			final List<String> original, final Case[] cases, final boolean flip, final boolean fixed) throws Exception{
		final List<String> expected=fixed ? splitSamRows(cases, flip) : original;
		Files.write(work.resolve("expected-sam-route.sam"), expected, StandardCharsets.US_ASCII);
		final Path output=work.resolve("sam-route.sam");
		final int exit=run(launcher, input, output, work, "sam-route", flip);
		check(exit==0, "SAM route exit="+exit+"; inspect "+work.resolve("sam-route.stderr"));
		check(Files.size(work.resolve("sam-route.stdout"))==0, "SAM route wrote unexpected stdout.");
		final ArrayList<String> actual=new ArrayList<String>();
		for(final String line : Files.readAllLines(output, StandardCharsets.US_ASCII)){
			if(!line.isEmpty() && line.charAt(0)!='@'){actual.add(line);}
		}
		check(actual.size()==expected.size(), "SAM route expected"+expected.size()+" records; got "+actual.size());
		for(int i=0; i<expected.size(); i++){
			check(expected.get(i).equals(actual.get(i)), "SAM route changed input row "+i+": "+actual.get(i));
		}
	}

	/** Four seven-base inputs become three fragments each; four three-base inputs keep their full SAM rows.
	 * The literal arrays describe stored input orientation; no production Read flags or output supply this oracle.
	 */
	private static List<String> splitSamRows(final Case[] cases, final boolean flip){
		final ArrayList<String> expected=new ArrayList<String>();
		for(final Case item : cases){
			if(item.shortRead){
				expected.add(item.inputLine());
			}else{
				final boolean reverse=item.reverse && (!item.mapped || !flip);
				final String[] sequences=reverse ? REVERSE_PARTS : FORWARD_PARTS;
				final String[] qualities=reverse ? REVERSE_QUALS : FORWARD_QUALS;
				for(int part=0; part<3; part++){
					expected.add(unmappedLine(item.name+"_"+(part+1), reverse ? 20 : 4,
							sequences[part], qualities[part]));
				}
			}
		}
		return expected;
	}

	/** The non-SAM output route splits the same inputs and exposes the stored Read orientation. */
	private static void fastqRoute(final Path launcher, final Path input, final Path work,
			final Case[] cases, final boolean flip) throws Exception{
		final StringBuilder expected=new StringBuilder();
		for(final Case item : cases){
			if(item.shortRead){appendFastq(expected, item.name, item.internalSequence(flip), item.internalQuality(flip));}
			else{
				final boolean reverse=item.reverseStored(flip);
				for(int i=0; i<3; i++){
					appendFastq(expected, item.name+"_"+(i+1), (reverse ? REVERSE_PARTS : FORWARD_PARTS)[i],
							(reverse ? REVERSE_QUALS : FORWARD_QUALS)[i]);
				}
			}
		}
		Files.write(work.resolve("expected-fastq-route.fq"), expected.toString().getBytes(StandardCharsets.US_ASCII));
		final Path output=work.resolve("fastq-route.fq");
		final int exit=run(launcher, input, output, work, "fastq-route", flip);
		check(exit==0, "FASTQ route exit="+exit+"; inspect "+work.resolve("fastq-route.stderr"));
		check(Files.size(work.resolve("fastq-route.stdout"))==0, "FASTQ route wrote unexpected stdout.");
		final String actual=new String(Files.readAllBytes(output), StandardCharsets.US_ASCII);
		check(expected.toString().equals(actual), "FASTQ route differs; inspect expected/actual files under "+work);
	}

	/** Runs the actual launcher under GNU timeout's separate process-group supervision.
	 * Normal exits pass through; timeout/forced-cleanup exits fail the route's exit check.
	 */
	private static int run(final Path launcher, final Path input, final Path output,
			final Path work, final String stem, final boolean flip) throws Exception{
		final List<String> command=Arrays.asList("/usr/bin/timeout", "--verbose", "--signal=TERM", "--kill-after=2s", "30s",
				"bash", launcher.toString(), "-Xmx128m", "-Xms32m",
				"in="+input, "out="+output, "t=1", "wt=1", "tin=1", "tout=1", "int=f", "ow=f",
				"simd=f", "minlength=1", "breaklength=3", "flipsam="+flip, "qin=33", "qout=33");
		Files.write(work.resolve(stem+".args"), command, StandardCharsets.UTF_8);
		final Process process=new ProcessBuilder(command).directory(launcher.getParent().toFile())
				.redirectOutput(work.resolve(stem+".stdout").toFile())
				.redirectError(work.resolve(stem+".stderr").toFile()).start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(60, TimeUnit.SECONDS)){
				Files.write(work.resolve(stem+".timeout"), Collections.singletonList("Owned launcher exceeded60 seconds."), StandardCharsets.US_ASCII);
				throw new AssertionError("Finite control exceeded60 seconds; inspect "+work.resolve(stem+".stderr"));
			}
			final int exit=process.exitValue();
			Files.write(work.resolve(stem+".exit"), Collections.singletonList(Integer.toString(exit)), StandardCharsets.US_ASCII);
			return exit;
		}catch(Exception e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
		}
	}

	/** Stops and reaps the exact owned GNU timeout supervisor.
	 * Its process group covers the launcher and Java child; allow its two-second KILL grace first.
	 */
	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 5)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){throw new AssertionError("Owned launcher did not terminate after forced cleanup.");}
	}

	/** Bounds cleanup even during repeated interrupts, restoring interruption after the owned wait. */
	private static boolean waitOwnedChild(final Process process, final int seconds){
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(seconds);
		boolean interrupted=false;
		try{
			while(process.isAlive()){
				final long remaining=deadline-System.nanoTime();
				if(remaining<=0){return false;}
				try{return process.waitFor(remaining, TimeUnit.NANOSECONDS);}
				catch(InterruptedException e){interrupted=true;}
			}
			return true;
		}finally{if(interrupted){Thread.currentThread().interrupt();}}
	}

	/** Pins parsing and default unmapped-tag generation; SavedStatics restores these configuration fields. */
	private static void configureSam(final boolean flip){
		SamLine.FLIP_ON_LOAD=flip;
		SamLine.PARSE_0=SamLine.PARSE_2=SamLine.PARSE_5=SamLine.PARSE_6=true;
		SamLine.PARSE_7=SamLine.PARSE_8=SamLine.PARSE_10=SamLine.PARSE_OPTIONAL=true;
		SamLine.PARSE_OPTIONAL_MD_ONLY=SamLine.PARSE_OPTIONAL_MATEQ_ONLY=false;
		SamLine.RNAME_AS_BYTES=true;
		SamLine.MAKE_CUSTOM_TAGS=SamLine.MAKE_TIME_TAG=false;
		SamLine.READGROUP_ID=null;
		SamLine.NO_TAGS=false;
		ReadStreamWriter.USE_ATTACHED_SAMLINE=true;
		FASTQ.PARSE_CUSTOM=false;
	}

	/** Eight literal input cases: four mapping/strand combinations, each long and exactly three bases. */
	private static Case[] cases(){
		return new Case[]{new Case("mp_long", true, false, false), new Case("mr_long", true, true, false),
			new Case("up_long", false, false, false), new Case("ur_long", false, true, false),
			new Case("mp_short", true, false, true), new Case("mr_short", true, true, true),
			new Case("up_short", false, false, true), new Case("ur_short", false, true, true)};
	}

	private static String unmappedLine(final String name, final int flag, final String sequence, final String quality){
		return name+"\t"+flag+"\t*\t0\t0\t*\t*\t0\t0\t"+sequence+"\t"+quality;
	}
	private static void appendFastq(final StringBuilder out, final String name, final String sequence, final String quality){
		out.append('@').append(name).append('\n').append(sequence).append("\n+\n").append(quality).append('\n');
	}
	private static String bases(final Read read){return new String(read.bases, StandardCharsets.US_ASCII);}
	private static String qualities(final Read read){
		final StringBuilder out=new StringBuilder();
		for(final byte quality : read.quality){out.append((char)(quality+33));}
		return out.toString();
	}
	private static String details(final Read read){
		if(read==null){return "null";}
		return read.id+"\tbases="+bases(read)+"\tquality="+qualities(read)+"\tmapped="+read.mapped()
				+"\tstrand="+read.strand()+"\tchrom="+read.chrom+"\tstart="+read.start+"\tstop="+read.stop
				+"\tnumericID="+read.numericID+"\tattached="+(read.samline!=null);
	}
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Literal stored SAM bytes and independent baseline parsing expectations for one input. */
	private static final class Case{
		Case(final String name_, final boolean mapped_, final boolean reverse_, final boolean shortRead_){
			this(name_, mapped_, reverse_, shortRead_, null);
		}
		Case(final String name_, final boolean mapped_, final boolean reverse_, final boolean shortRead_, final Integer x6_){
			name=name_; mapped=mapped_; reverse=reverse_; shortRead=shortRead_;
			x6=x6_;
		}
		String inputLine(){
			final String sequence=reverse ? (shortRead ? "CGT" : "GTAACGT") : (shortRead ? "ACG" : "ACGTTAC");
			final String quality=reverse ? (shortRead ? "CBA" : "GFEDCBA") : (shortRead ? "ABC" : "ABCDEFG");
			return line(sequence, quality);
		}
		/** rcomp changes shared short arrays but leaves the attached SAM flag and other fields intact. */
		String shortOutput(final boolean rcomp){
			return rcomp ? line(reverse ? "ACG" : "CGT", reverse ? "ABC" : "CBA") : inputLine();
		}
		private String line(final String sequence, final String quality){
			final String tags="\tZZ:Z:keep_"+name+(x6==null ? "" : "\tX6:i:"+x6);
			if(!mapped){return unmappedLine(name, reverse ? 20 : 4, sequence, quality)+tags;}
			return name+"\t"+(reverse ? 16 : 0)+"\tchr1\t"+(shortRead ? 21 : 11)
					+"\t"+(shortRead ? 37 : 42)+"\t"+(shortRead ? "3M" : "7M")
					+"\t*\t0\t0\t"+sequence+"\t"+quality
					+"\tNM:i:0\tMD:Z:"+(shortRead ? 3 : 7)+"\tAS:i:"+(shortRead ? 3 : 7)+tags;
		}
		boolean reverseStored(final boolean flip){return reverse && (!mapped || !flip);}
		//X6 uses internal Read flags: low bit is strand, next bit is mapped. Keep the oracle literal.
		boolean readReverse(){return x6==null ? reverse : (x6.intValue()&1)!=0;}
		boolean readMapped(){return x6==null ? mapped : (x6.intValue()&2)!=0;}
		String internalSequence(final boolean flip){return reverseStored(flip) ? (shortRead ? "CGT" : "GTAACGT") : (shortRead ? "ACG" : "ACGTTAC");}
		String internalQuality(final boolean flip){return reverseStored(flip) ? (shortRead ? "CBA" : "GFEDCBA") : (shortRead ? "ABC" : "ABCDEFG");}
		final String name;
		final boolean mapped, reverse, shortRead;
		final Integer x6;
	}

	/** Expected serialized text and original object identities, built before production splitting. */
	private static final class Expected{
		Expected(final String row_, final Read original_, final SamLine line_, final long numericID_, final boolean split_,
				final boolean mapped_, final boolean reverse_){
			row=row_; original=original_; line=line_; numericID=numericID_; split=split_;
			mapped=mapped_; reverse=reverse_;
		}
		final String row;
		final Read original;
		final SamLine line;
		final long numericID;
		final boolean split, mapped, reverse;
	}

	/** Records one synchronous driver submission; production SamWriter supplies the actual serializer. */
	private static final class RecordingWriter implements Writer{
		public void add(final ArrayList<Read> list, final long id){addReads(new ListNum<Read>(list, id));}
		public void addReads(final ListNum<Read> batch){ids.add(batch.id); reads.addAll(batch.list);}
		public void addLines(final ListNum<SamLine> lines){throw new UnsupportedOperationException("Read-route capture only.");}
		public void start(){}
		public long readsWritten(){return reads.size();}
		public long basesWritten(){long count=0; for(final Read read : reads){if(read!=null){count+=read.length();}} return count;}
		public void poison(){}
		public boolean waitForFinish(){return false;}
		public boolean poisonAndWait(){return false;}
		public void finishError(){throw new AssertionError("Unexpected writer abandonment.");}
		public String fname(){return "recording-read-route";}
		public boolean errorState(){return false;}
		public boolean finishedSuccessfully(){return true;}
		final ArrayList<Long> ids=new ArrayList<Long>();
		final ArrayList<Read> reads=new ArrayList<Read>();
	}

	/** Restores mutable configuration fields in the known owners after constructor/parser use.
	 * No production workers are started in the direct path. Fresh JVM isolation also
	 * contains initialization caches and the processor's thread-local random state.
	 */
	private static final class SavedStatics{
		SavedStatics(final Class<?>... types) throws IllegalAccessException{
			for(final Class<?> type : types){
				for(final Field field : type.getDeclaredFields()){
					final int flags=field.getModifiers();
					if(Modifier.isStatic(flags) && !Modifier.isFinal(flags)){
						field.setAccessible(true);
						fields.add(field); values.add(field.get(null));
					}
				}
			}
		}
		void restore() throws IllegalAccessException{
			for(int i=fields.size()-1; i>=0; i--){fields.get(i).set(null, values.get(i));}
		}
		final ArrayList<Field> fields=new ArrayList<Field>();
		final ArrayList<Object> values=new ArrayList<Object>();
	}

	private static final String[] FORWARD_PARTS={"ACG", "TTA", "C"}, FORWARD_QUALS={"ABC", "DEF", "G"};
	private static final String[] REVERSE_PARTS={"GTA", "ACG", "T"}, REVERSE_QUALS={"GFE", "DCB", "A"};
	private static final String[] RC_FORWARD_PARTS={"CGT", "TAA", "G"}, RC_FORWARD_QUALS={"CBA", "FED", "G"};
	private static final String[] RC_REVERSE_PARTS={"TAC", "CGT", "A"}, RC_REVERSE_QUALS={"EFG", "BCD", "A"};
	private static int checks;
}
