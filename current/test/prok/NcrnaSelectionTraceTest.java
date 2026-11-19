package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReportedOrf;
import prok.NcrnaSelectionTrace;
import prok.NcrnaSelectionAudit;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import prok.ProkObject;
import shared.Shared;
import structures.ByteBuilder;

/** Independent literal overlaps on selected edges and newly overlapping nonedges.
 * The observer must preserve every path object and coordinate on both strands.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSelectionTraceTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaSelectionTraceTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown selected-geometry fixture argument: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1 && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			for(int strand=0; strand<2; strand++){run(new File(parser.out1), strand);}
			twoTargets(new File(parser.out1));
			physicalMovements(new File(parser.out1));
			pp.outstream.println("NCRNA_SELECTION_TRACE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Every selected call contributes, including unchanged third-member calls. */
	private static void physicalMovements(final File root){
		final String before=new File(root, "movement_before.tsv").getPath();
		final String after=new File(root, "movement_after.tsv").getPath();
		final String[] members={"cobalamin_rf00174", "cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
		for(int mode=0; mode<2; mode++){
			try(Output output=new Output(mode==0 ? before : after)){
				output.writer.print(NcrnaSelectionTrace.HEADER);
				output.writer.println("I\t0\tliteral_movement.fa");
				for(int strand=0; strand<2; strand++){
					final byte[] bases=new byte[1000];
					final ArrayList<Orf> calls=new ArrayList<Orf>();
					for(int i=0; i<4; i++){
						final int start=100+150*i;
						final int stop=start+59;
						final NcrnaReportedOrf call=new NcrnaReportedOrf("literal", start, stop, strand, bases);
						call.ncrnaFamily=members[i];
						call.trnaModel=members[i]+"_literal";
						call.orfScore=42;
						if(mode==1){
							if(i==0 || i==2){call.start-=2;}
							if(i==1 || i==2){call.stop+=3;}
							call.captureRefinement(start, stop, 16, 1000);
						}
						if(strand==1){call.flip();}
						calls.add(call);
					}
					NcrnaSelectionTrace.write(output.writer, new ByteBuilder(), 0, strand, "strand"+strand,
						"cobalamin", calls, 80, 110);
				}
			}
		}
		final String summary=new File(root, "movement_summary.tsv").getPath();
		final String movement=new File(root, "movement_counts.tsv").getPath();
		NcrnaSelectionAudit.audit(before, after, summary, movement, 1);
		final String[] expected={
			"before\tcobalamin_rf00174\t4\t0\t0\t0\t0",
			"before\tcobalamin_rf01482\t2\t0\t0\t0\t0",
			"before\tcobalamin_rf01689\t2\t0\t0\t0\t0",
			"before\tcobalamin\t8\t0\t0\t0\t0",
			"after\tcobalamin_rf00174\t4\t2\t2\t0\t4",
			"after\tcobalamin_rf01482\t2\t2\t2\t2\t2",
			"after\tcobalamin_rf01689\t2\t0\t0\t0\t0",
			"after\tcobalamin\t8\t4\t4\t2\t6"};
		try(Input input=new Input(movement)){
			check(new String(input.reader.nextLine(), java.nio.charset.StandardCharsets.US_ASCII)
				.equals("view\tmember\tselected_calls\tmoved5\tmoved3\tmoved_both\tmoved_either"), "Movement schema must state the physical-call denominator");
			for(final String row:expected){
				check(java.util.Arrays.equals(row.getBytes(java.nio.charset.StandardCharsets.US_ASCII), input.reader.nextLine()),
					"Literal plus/minus, one/both/unchanged movement counts differ: "+row);
			}
			check(input.reader.nextLine()==null, "Exactly four cohorts in each view are required");
		}
		final String historical=new File(root, "movement_historical_summary.tsv").getPath();
		NcrnaSelectionAudit.audit(before, after, historical, 1);
		try(Input left=new Input(summary); Input right=new Input(historical)){
			for(byte[] row; (row=left.reader.nextLine())!=null;){
				check(java.util.Arrays.equals(row, right.reader.nextLine()), "Optional movement reporting must preserve every old summary byte");
			}
			check(right.reader.nextLine()==null, "Historical summary cannot gain or lose rows");
		}
	}

	/** Two moved RNAs share one pair:20 ordinary overlap bases become41, crossing30. */
	private static void twoTargets(final File root){
		final String off=new File(root, "cap_off.tsv").getPath(), on=new File(root, "cap_on.tsv").getPath();
		for(int mode=0; mode<2; mode++){
			final byte[] bases=new byte[240];
			final NcrnaReportedOrf first=new NcrnaReportedOrf("literal", 100, 159, 0, bases);
			final NcrnaReportedOrf second=new NcrnaReportedOrf("literal", 140, 199, 0, bases);
			first.ncrnaFamily=second.ncrnaFamily="cobalamin_rf00174";
			first.trnaModel=second.trnaModel="cobalamin_rf00174_literal";
			if(mode==1){
				first.stop=170;
				first.captureRefinement(100, 159, 16, 1000);
				second.start=130;
				second.captureRefinement(140, 199, 16, 1000);
			}
			second.prevPlus=first;
			final ArrayList<Orf> path=new ArrayList<Orf>();
			path.add(first);
			path.add(second);
			try(Output output=new Output(mode==0 ? off : on)){
				output.writer.print(NcrnaSelectionTrace.HEADER);
				output.writer.println("I\t0\tliteral.fa");
				if(mode==1){NcrnaSelectionTrace.write(output.writer, new ByteBuilder(), 0, 1, "empty", "cobalamin", new ArrayList<Orf>(), 30, 40);}
				NcrnaSelectionTrace.write(output.writer, new ByteBuilder(), 0, 0, "literal", "cobalamin", path, 30, 40);
				if(mode==0){NcrnaSelectionTrace.write(output.writer, new ByteBuilder(), 0, 1, "empty", "cobalamin", new ArrayList<Orf>(), 30, 40);}
			}
		}
		final LineParser1 p=new LineParser1('\t');
		int pairs=0;
		try(Input input=new Input(on)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(row.length<2 || row[0]!='Q'){continue;}
				p.set(row);
				check(p.parseInt(4)==0 && p.parseInt(5)==1 && p.parseInt(9)==20 && p.parseInt(10)==41
					&& p.parseInt(11)==30 && p.parseInt(6)==1, "Two moved targets must emit their one literal selected edge exactly once");
				pairs++;
			}
		}
		check(pairs==1, "The reverse target iteration must not emit the same pair again");
		final String summary=new File(root, "cap_summary.tsv").getPath();
		NcrnaSelectionAudit.audit(off, on, summary, 1);
		metric(summary, "exact_path_multiset", 1);
		metric(summary, "after_all_pairs_increased_bases", 21);
		metric(summary, "after_all_pairs_cap_crossings", 1);
		metric(summary, "after_all_pairs_affected_target_rnas", 2);
		final String changed=new File(root, "changed_score.tsv").getPath();
		final String duplicate=new File(root, "duplicate_pair.tsv").getPath();
		mutate(on, changed, false);
		mutate(on, duplicate, true);
		fails(off, changed, new File(root, "changed_summary.tsv").getPath());
		fails(off, duplicate, new File(root, "duplicate_summary.tsv").getPath());
		for(int corruption=0; corruption<3; corruption++){
			final String bad=new File(root, "corrupt_"+corruption+".tsv").getPath();
			corrupt(on, bad, corruption);
			final String result=new File(root, "corrupt_summary_"+corruption+".tsv").getPath();
			fails(off, bad, result);
			if(corruption==2){metric(result, "exact_path_multiset", 0);}
		}
	}

	/** Deleted pair, in-range wrong order, and an extra identical P identity are distinct failures. */
	private static void corrupt(final String in, final String out, final int mode){
		final LineParser1 p=new LineParser1('\t');
		boolean changed=false;
		try(Input input=new Input(in); Output output=new Output(out)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(row.length<2 || row[0]=='#'){
					output.writer.println(row);
					continue;
				}
				p.set(row);
				if(row[0]=='Q' && mode<2){
					if(mode==1){field(output, p, row, 12, 0);}
					changed=true;
				}else if(row[0]=='P' && mode==2 && p.parseInt(4)==1){
					output.writer.println(row);
					field(output, p, row, 4, 2);
					changed=true;
				}else if(row[0]=='T' && mode==2 && p.parseInt(4)==2){field(output, p, row, 4, 3);}
				else{output.writer.println(row);}
			}
		}
		check(changed, "The corrupt-evidence fixture must actually alter the intended record");
	}
	private static void field(final Output output, final LineParser1 p, final byte[] row, final int field, final int value){
		p.setBounds(field);
		output.writer.print(new ByteBuilder().append(row, 0, p.a()).append(value).append(row, p.b(), row.length-p.b()).nl());
	}

	/** An independently changed score is a path mismatch; a repeated Q row is corrupt evidence. */
	private static void mutate(final String in, final String out, final boolean duplicate){
		boolean changed=false;
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		try(Input input=new Input(in); Output output=new Output(out)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(!duplicate && !changed && row.length>1 && row[0]=='P'){
					p.set(row);
					p.setBounds(14);
					output.writer.print(b.clear().append(row, 0, p.a()).append(1).nl());
					changed=true;
				}else{
					output.writer.println(row);
					if(duplicate && row.length>1 && row[0]=='Q'){
						output.writer.println(row);
						changed=true;
					}
				}
			}
		}
		check(changed, "The negative must actually change the intended native record");
	}
	private static void metric(final String path, final String key, final long expected){
		final LineParser1 p=new LineParser1('\t');
		boolean found=false;
		try(Input input=new Input(path)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals(key, 0)){
					check(!found && p.parseLong(1)==expected, "Wrong literal selected-geometry summary for "+key);
					found=true;
				}
			}
		}
		check(found, "Missing literal selected-geometry metric: "+key);
	}
	private static void fails(final String before, final String after, final String out){
		try{NcrnaSelectionAudit.audit(before, after, out, 1);}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Changed path/corrupt pair evidence must fail");
	}

	private static void run(final File root, final int strand){
		final byte[] bases=new byte[240];
		final NcrnaReportedOrf rna=new NcrnaReportedOrf("literal", 100, 159, strand, bases);
		rna.ncrnaFamily="cobalamin_rf00174";
		rna.trnaModel="cobalamin_rf00174_literal";
		rna.orfScore=42;
		rna.stop=162;
		rna.captureRefinement(100, 159, 16, 1000);
		if(strand==1){rna.flip();}
		final Orf cds=new Orf("literal", strand==0 ? 120 : 60, strand==0 ? 220 : 100,
			strand, 1, bases, false, ProkObject.CDS);
		final Orf other=new Orf("literal", strand==0 ? 161 : 78, strand==0 ? 161 : 78,
			1-strand, 0, bases, false, ProkObject.RNA);
		other.ncrnaFamily="sam1";
		other.trnaModel="sam1_literal";
		if(strand==0){cds.prevPlus=rna;}
		else{rna.prevPlus=cds;}
		final ArrayList<Orf> path=new ArrayList<Orf>();
		path.add(rna);
		path.add(cds);
		path.add(other);
		final int pathStart=rna.start, pathStop=rna.stop;
		final String output=new File(root, "geometry_"+strand+".tsv").getPath();
		try(Output out=new Output(output)){
			out.writer.print(NcrnaSelectionTrace.HEADER);
			NcrnaSelectionTrace.write(out.writer, new ByteBuilder(), 0, 0, "literal", "cobalamin", path, 80, 110);
			NcrnaSelectionTrace.write(out.writer, new ByteBuilder(), 0, 1, "empty", "cobalamin", new ArrayList<Orf>(), 80, 110);
		}
		check(path.size()==3 && path.get(0)==rna && rna.start==pathStart && rna.stop==pathStop && rna.orfScore==42,
			"Observation must leave membership, path coordinates and cached score untouched");
		final LineParser1 p=new LineParser1('\t');
		int selected=0, pairs=0, edges=0, newlyOverlapping=0, increase=0, orderChanges=0, totals=0;
		try(Input input=new Input(output)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				if(p.termEquals("P", 0)){
					check(p.terms()==15, "Selected records preserve the complete declared geometry schema");
					selected++;
				}else if(p.termEquals("Q", 0)){
					check(p.terms()==16 && p.parseInt(4)==0 && p.parseInt(5)>0, "Pair rows must identify the moved RNA and its exact neighbour once");
					final int before=p.parseInt(9), after=p.parseInt(10);
					if(p.parseInt(5)==1){
						check(before==(strand==0 ? 40 : 21) && after==(strand==0 ? 43 : 24)
							&& p.parseInt(6)==1 && p.parseInt(11)==80,
							"Actual selected CDS edge gains exactly three overlap bases below its unchanged cap");
					}else{
						check(before==0 && after==1 && p.parseInt(6)==0 && p.parseInt(11)==110,
							"A separate nonedge gains one overlap base on the opposite strand");
					}
					edges+=p.parseInt(6);
					newlyOverlapping+=before==0 && after>0 ? 1 : 0;
					increase+=after-before;
					orderChanges+=p.parseInt(12)==p.parseInt(13) ? 0 : 1;
					pairs++;
				}else if(p.termEquals("T", 0)){
					check(p.terms()==6 && p.parseInt(4)==(totals==0 ? 3 : 0) && p.parseInt(5)==(totals==0 ? 2 : 0),
						"Empty contigs retain explicit zero selected-feature and pair counts");
					totals++;
				}else{throw new AssertionError("Unexpected selected-geometry record");}
			}
		}
		check(selected==3 && pairs==2 && edges==1 && newlyOverlapping==1 && increase==4 && totals==2,
			"Literal selected membership, evaluated edge and nonedge overlap accounting must conserve both geometries");
		check(orderChanges==(strand==0 ? 1 : 0), "Only the forward singleton crosses stop order in this literal fixture");
	}

	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
}
