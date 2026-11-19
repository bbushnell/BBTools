package test.prok;

import java.io.File;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Arrays;
import java.nio.charset.StandardCharsets;

import map.LongHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.MitoSsuShortlistControl;
import prok.NcrnaStageDiagLogSink;
import prok.NcrnaStageDiagSink;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal presence and framing oracles; optional old/new real diagnostic comparison.
 * Tiny finite alignments check producer presence; no fitting or old-zero inference.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaStageDiagTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaStageDiagTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String old=null;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
				if(a.equals("old")){old=b;}
				else{check(parser.parse(arg, a, b), "Use old=baseline in=current out=fresh_report, or only out=fresh_fixture_directory");}
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && (old==null)==(parser.in1==null), "Bind a fresh fixture or both diagnostic versions");
			if(old==null){unit(parser.out1);}else{compare(old, parser.in1, parser.out1);}
			pp.outstream.println("NCRNA_STAGE_DIAGNOSTIC_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void unit(final String directory){
		final File root=new File(directory);
		check(root.mkdir(), "Stage fixture directory must be fresh");
		final String current=new File(root, "current.tsv").getPath(), old=new File(root, "old.tsv").getPath();
		try(Output out=new Output(current)){
			out.writer.println(NcrnaStageDiagLogSink.HEADER);
			final NcrnaStageDiagLogSink sink=new NcrnaStageDiagLogSink(out.writer);
			sink.seed("16S", "contig", 0, 100, new int[0]);
			sink.seed("16S", "contig", 0, 100, new int[]{8, 20});
			sink.model("16S", "contig", 0, 100, 1, 0, 99, 0, "-1", 0);
			sink.window("16S", "contig", 0, 100, 1, 0, 99, "REJECT_ID", 0, 0, -1, null,
				0, 0, -999, -1, -1, 0, false);
			sink.window("16S", "contig", 0, 100, 1, 0, 99, "REJECT_ID", 0, 0, 0, "-1",
				0, 0, 0, 0, 99, NcrnaStageDiagSink.PRIMARY_ID|NcrnaStageDiagSink.RESCUE_ID|NcrnaStageDiagSink.HBM|NcrnaStageDiagSink.TRIM, false);
			sink.path("tmrna", "contig", 0, 100, "-1", 0, 99, false);
			sink.path("16S", "contig", 0, 100, "-1", 0, 99, true, 0);
			sink.model("16S", "contig", 0, 100, 1, 0, 99, 0, "tiny", .0001f);
			fails(()->sink.window("16S", "contig", 0, 100, 1, 0, 99, "REJECT_ID", 0, 0, 0, "-1",
				Float.NaN, 0, 0, 0, 99, NcrnaStageDiagSink.PRIMARY_ID, false));
		}
		final LineParser1 p=new LineParser1('\t');
		int row=0;
		try(Input in=new Input(current)){
			NcrnaStageDiagLogSink.readHeader(in.reader, false, current);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				p.set(line);
				NcrnaStageDiagLogSink.validateRow(p, false);
				check(p.terms()==21, "A trailing empty seed-position field must survive native parsing");
				if(row==0){check(p.parseInt(19)==0 && p.length(20)==0 && p.termEquals("-1", 13), "Zero observed seeds is real; unavailable identity and positions are distinct");}
				if(row==1){check(p.parseInt(19)==2 && p.termEquals("8,20", 20), "Observed positions retain their literal order");}
				if(row==2){check(p.termEquals("-1", 12) && p.parseFloat(13)==0 && p.termEquals("-1", 14), "A model named-1 and its measured zero identity remain present");}
				if(row==3){check(p.termEquals("-1", 13) && p.termEquals("-1", 14) && p.termEquals("-1", 15) && p.length(18)==0, "Unmeasured scores and unattempted trim must not become zero/false");}
				if(row==4){check(p.parseFloat(13)==0 && p.parseFloat(14)==0 && p.parseFloat(15)==0 && p.termEquals("false", 18), "Every measured zero and actual failed trim is retained");}
				if(row==5){check(p.termEquals("-1", 13) && p.termEquals("-1", 12), "A missing generic verifier identity does not erase a literal model name-1");}
				if(row==6){check(p.parseFloat(13)==0 && p.termEquals("DP_RETAINED", 8), "Actual accepted verifier identity zero survives PATH output");}
				if(row==7){check(NcrnaStageDiagLogSink.score(p, 13)==.0001f, "Scientific notation must retain an actual small measurement");}
				row++;
			}
		}
		check(row==8, "A failed nonfinite write must not emit a partial event");
		try(Output out=new Output(old)){
			out.writer.println(NcrnaStageDiagLogSink.LEGACY_HEADER);
			out.writer.println("MODEL\t16S\tcontig\t0\t100\t1\t0\t99\tALIGNMENT\t.\t.\t0\t-1\t0.0\t.\t.\t.\t.\t.\t.\t.");
		}
		fails(()->{try(Input in=new Input(old)){NcrnaStageDiagLogSink.readHeader(in.reader, false, old);}});
		check(MitoSsuShortlistControl.rows(old, "contig", 0, true).size()==1, "Historical rows are read only after explicit selection");
		final String prefix=new File(root, "must_not_exist").getPath();
		fails(()->MitoSsuShortlistControl.main(new String[]{"stages="+old, "cases=/not_a_resource", "outprefix="+prefix, "historicalstages=t"}));
		check(!new File(prefix+".top1.stages.tsv").exists(), "Historical raw-row parity must fail before resource loading/output/caller execution");
		producerPrimary(root, false);
		producerPrimary(root, true);
		integerFields();
	}
	/** Literal overflow must not wrap into measured zero or the absence marker. */
	private static void integerFields(){
		final String prefix="SEED\ttest\tcontig\t0\t100\t-1\t-1\t-1\t\t-1\t-1\t-1\t\t-1\t-1\t-1\t-1\t-1\t\t";
		final LineParser1 parser=new LineParser1('\t');
		parser.set((prefix+"0\t").getBytes(StandardCharsets.US_ASCII));
		NcrnaStageDiagLogSink.validateRow(parser, false);
		for(final String malformed:new String[]{"4294967296", "4294967295", "2147483648", "", "2x", "-01"}){
			parser.set((prefix+malformed+'\t').getBytes(StandardCharsets.US_ASCII));
			fails(()->NcrnaStageDiagLogSink.validateRow(parser, false));
		}
	}

	/** Scrabble measures zero without a winner. Quantum requires a nonempty
	 * traceback, so its finite fixture uses an exact match with real identity one. */
	private static void producerPrimary(final File root, final boolean reused){
		final byte[] model=new byte[64], genome=new byte[64];
		Arrays.fill(model, (byte)'A');
		Arrays.fill(genome, reused ? (byte)'A' : (byte)'C');
		final LongHashSet seeds=new LongHashSet(1);
		seeds.add(reused ? 0 : 0x155555555L);//Actual forward A17 or C17.
		final NcrnaScavenger caller=new NcrnaScavenger(new byte[][]{model}, null, new String[]{"zero"},
			seeds, 17, 20, 64, 7, 1, false, 0, 0, 0, 1, 0, 1, .9f, .8f);
		caller.setMappingOracleExhaustive(true);
		caller.scavengePass2=false;
		caller.reuseConsensusAlignment=reused;
		caller.trimAlignmentExtent=false;
		final String path=new File(root, reused ? "producer_reused.tsv" : "producer_separate.tsv").getPath();
		try(Output out=new Output(path)){
			out.writer.println(NcrnaStageDiagLogSink.HEADER);
			caller.setStageDiagSink(new NcrnaStageDiagLogSink(out.writer), "test");
			check(caller.scavenge("contig", genome, 0, new ArrayList<int[]>(), new int[]{20}).size()==(reused ? 1 : 0),
				"Exact Quantum match is accepted; unrelated Scrabble sequences are rejected");
		}
		final LineParser1 parser=new LineParser1('\t');
		int windows=0, models=0;
		try(Input in=new Input(path)){
			NcrnaStageDiagLogSink.readHeader(in.reader, false, path);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				NcrnaStageDiagLogSink.validateRow(parser, false);
				if(parser.termEquals("WINDOW", 0)){
					check(parser.parseFloat(13)==(reused ? 1 : 0) && parser.parseInt(11)==(reused ? 0 : -1) && parser.length(18)==0,
						"Measured identity survives with or without a winner; unattempted trim stays empty");
					windows++;
				}else if(parser.termEquals("MODEL", 0)){
					check(parser.parseFloat(13)==(reused ? 1 : 0), "Actual model identities must retain exact-match one and mismatch zero");
					models++;
				}
			}
		}
		check(windows==1 && models==1, "Both primary implementations must actually align and emit their measured identity");
	}

	/** Current absence is accepted only against the old field's declared absence
	 * spelling. All present values and event order remain byte-identical. */
	private static void compare(final String old, final String current, final String output){
		check(!new File(output).exists(), "Comparison report must be fresh");
		final LineParser1 a=new LineParser1('\t'), b=new LineParser1('\t');
		long rows=0, unchanged=0, absent=0;
		final long[] kinds=new long[4];
		try(Input before=new Input(old); Input after=new Input(current)){
			NcrnaStageDiagLogSink.readHeader(before.reader, true, old);
			NcrnaStageDiagLogSink.readHeader(after.reader, false, current);
			for(byte[] first; (first=before.reader.nextLine())!=null;){
				final byte[] second=after.reader.nextLine();
				check(second!=null, "The new diagnostic cannot lose events");
				a.set(first);
				b.set(second);
				NcrnaStageDiagLogSink.validateRow(a, true);
				NcrnaStageDiagLogSink.validateRow(b, false);
				for(int column=0; column<21; column++){
					if(sameField(a, b, column)){
						unchanged++;
						continue;
					}
					final boolean numeric=b.termEquals("-1", column) && (a.termEquals(".", column)
						|| column>=13 && column<=15 && a.parseFloat(column)==0
						|| column==15 && a.termEquals("-999.0", column));
					final boolean text=b.length(column)==0 && (a.termEquals(".", column) || column==18 && a.termEquals("false", column));
					check(numeric || text, "A measured value/event changed: row="+rows+" column="+column);
					absent++;
				}
				final int kind=b.termEquals("SEED", 0) ? 0 : b.termEquals("WINDOW", 0) ? 1 : b.termEquals("MODEL", 0) ? 2 : 3;
				kinds[kind]++;
				rows++;
			}
			check(after.reader.nextLine()==null && absent>0, "The comparison must conserve events and exercise absent-slot changes");
		}
		for(final long count:kinds){check(count>0, "The real comparison must cover all four event types");}
		try(Output out=new Output(output)){
			out.writer.println("events\tunchanged_fields\tchanged_absent_fields\tseed\twindow\tmodel\tpath");
			final ByteBuilder row=new ByteBuilder().append(rows).tab().append(unchanged).tab().append(absent);
			for(final long count:kinds){row.tab().append(count);}
			out.writer.print(row.nl());
		}
	}
	/** Compare retained byte ranges without allocating a String for each field. */
	private static boolean sameField(final LineParser1 first, final LineParser1 second, final int column){
		final int length=first.setBounds(column);
		if(second.setBounds(column)!=length){return false;}
		final byte[] a=first.line(), b=second.line();
		for(int i=0; i<length; i++){
			if(a[first.a()+i]!=b[second.a()+i]){return false;}
		}
		return true;
	}
	private static void fails(final Runnable action){
		boolean failed=false;
		try{action.run();}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "The malformed or unselected historical format must fail loudly");
	}
	private static void check(final boolean valid, final String message){if(!valid){throw new AssertionError(message);}}
}
