package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaGateCutoffSummary;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Synthetic saved tables distinguish medians, objective ties and outer leakage.
 * @author Brian, Ganyu, Raiden
 */
public final class NcrnaGateCutoffSummaryTest{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaGateCutoffSummaryTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Use native flag=value fixture options");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append,
				"Use out=fresh_directory");
			final File root=new File(parser.out1);
			check(root.mkdir(), "Fixture directory must be fresh");
			run(root);
			pp.outstream.println("NCRNA_GATE_CUTOFF_SUMMARY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}
	/*--------------------------------------------------------------*/
	/*----------------          Finite Cases        ----------------*/
	/*--------------------------------------------------------------*/

	/** Full synthetic groups expose a wrong mean, lost ties or outer-selected cutoff. */
	private static void run(final File root){
		final File outer=new File(root, "outer");
		check(outer.mkdir(), "Fresh outer fixture");
		write(new File(outer, "complete.txt"), "THIRTY_SIX_OUTER_GATE_COMPARISONS_COMPLETE_NO_SELECTION_NO_CALLER\n");
		final ByteBuilder manifest=new ByteBuilder("cell\tarm\tend\tseed\toutput\trate\trole\tjob\n");
		String first=null;
		for(final String end:new String[]{"5prime", "3prime"}){
			for(final String activation:new String[]{"sigmoid", "rslog"}){
				for(final String rate:new String[]{"0.003", "0.0003"}){
					for(int seed=58; seed<=60; seed++){
						final String cell="R_"+end+"_"+activation+"_lr"+rate+"_s"+seed;
						if(first==null){first=cell;}
						manifest.append(cell).append("\tR\t").append(end).tab().append(seed).tab().append(activation).tab().append(rate).append("\tROLE\tunused\n");
						for(int cohort=0; cohort<2; cohort++){
							final File run=new File(cohort==0 ? root : outer, "runs/"+cell);
							check(new File(run, "replay").mkdirs(), "Fresh synthetic model output");
							write(new File(run, "exit.txt"), "0\n");
							write(new File(run, "complete.txt"), cohort==0 ? "GATE_CONVERGENCE_FIT_AND_SAVED_REPLAY_COMPLETE_OUTER_UNPERFORMED\n"
								: "OUTER_GATE_SCORED_AND_REPLAYED_NO_SELECTION_NO_CALLER\n");
							final String prefix=cohort==0 ? "" : "fit.";
							for(final String name:new String[]{"source.txt", "build.seal", "job.txt"}){write(new File(run, prefix+name), "synthetic\n");}
							if(cohort==1){write(new File(run, "selection.txt"), "monitor_mse\n");}
							write(new File(run, "replay/summary.tsv"), table(end, seed, cohort==1));
						}
					}
				}
			}
		}
		final File cells=new File(root, "cells.tsv");
		write(cells, manifest.toString());
		collect(root, "good");
		final String choices=read(new File(root, "good/choices.tsv"));
		check(choices.contains("5prime\tsigmoid\t0.003\tmonitor\tmax_exact\tgate\t0.4\t20.0\t18.0\t2.0\t0.0\t120.0\t"),
			"Choose the median120, not the mean122, of asymmetric three-seed exactness");
		check(choices.contains("\tmonitor\tmax_exact\tgate\t0.7\t"), "Retain all monitor exactness ties");
		check(choices.contains("\tmonitor\tmax_net_improvements\tgate\t0.5\t"), "Net improvements use closer minus worse");
		check(choices.contains("\tmonitor\tmin_rms\tgate\t0.6\t"), "Minimum RMS is a separate objective");
		check(choices.contains("\touter\tmax_exact\tgate\t0.4\t") && !choices.contains("\tgate\t0.9\t"),
			"Outer's deliberately superior0.9 cutoff must never affect monitor selection");
		check(read(new File(root, "good/monitor_grid.tsv")).contains("\t1\t1\t1\t0\n"),
			"A nontrivial synthetic gate can beat fixed-rule extrema in exactness,RMS and net improvements");
		check(choices.contains("\tfixed_rule\tdo_nothing\t\t") && !choices.contains("\tNA\t"),
			"A fixed rule has no invented cutoff in normalized output");
		for(final String line:manifest.toString().split("\n")){
			if(line.startsWith("cell\t")){continue;}
			final String cell=line.substring(0, line.indexOf('\t'));
			for(final File base:new File[]{root, outer}){
				final File file=new File(base, "runs/"+cell+"/replay/summary.tsv");
				write(file, read(file).replace("\tNA\t", "\t\t"));
			}
		}
		collect(root, "current");
		for(final String name:new String[]{"monitor_grid.tsv", "choices.tsv", "summary.tsv"}){
			check(read(new File(root, "good/"+name)).equals(read(new File(root, "current/"+name))),
				"Historical and current missing-cutoff inputs must yield identical complete reports: "+name);
		}
		final File firstTable=new File(root, "runs/"+first+"/replay/summary.tsv");
		final String current=read(firstTable);
		write(firstTable, current.replace("do_nothing\t\t", "do_nothing\t0\t"));
		fails(()->collect(root, "invented_cutoff"));
		check(!new File(root, "invented_cutoff/summary.tsv").exists(), "A numeric fixed-rule cutoff must fail before completion");
		write(firstTable, current);
		write(new File(outer, "runs/"+first+"/fit.job.txt"), "different\n");
		fails(()->collect(root, "wrong_model"));
		write(new File(outer, "runs/"+first+"/fit.job.txt"), "synthetic\n");
		final String original=manifest.toString();
		write(cells, original.substring(0, original.lastIndexOf('\n', original.length()-2)+1));
		fails(()->collect(root, "missing_seed"));
		write(cells, original);
		final File table=new File(root, "runs/"+first+"/replay/summary.tsv");
		write(table, table("5prime", 58, false).replaceFirst("gate\t0.5[^\n]*\n", ""));
		fails(()->collect(root, "missing_cutoff"));
		for(final String name:new String[]{"wrong_model", "missing_seed", "missing_cutoff"}){
			check(!new File(root, name+"/summary.tsv").exists(), "Invalid evidence cannot publish completion");
		}
	}
	/** Distinct exactness/net/RMS optima, plus a deliberately tempting outer-only optimum. */
	private static String table(final String end, final int seed, final boolean outer){
		final ByteBuilder text=new ByteBuilder(NcrnaGateCutoffSummary.TABLE_HEADER).nl();
		for(int row=0; row<26; row++){
			final String[] values=new String[40];
			Arrays.fill(values, "0");
			final int cutoff=row<4 ? -1 : row<15 ? row-4 : row-15;
			values[0]=row<4 ? new String[]{"do_nothing", "margin0", "margin0.1", "all_proposals"}[row] : row<15 ? "gate" : "gate_hard_oor_veto";
			values[1]=cutoff<0 ? "NA" : Double.toString(cutoff/10.0);
			values[2]=end;
			final int trials=outer ? (end.equals("5prime") ? 13401 : 13406) : (end.equals("5prime") ? 88966 : 88979);
			values[3]=""+trials;
			values[4]=""+(outer ? (end.equals("5prime") ? 7616 : 9347) : (end.equals("5prime") ? 50183 : 61986));
			int moved=20, better=15, worse=4, exact=111;
			double rms=9;
			if(row==0){
				moved=better=worse=0;
				exact=100;
				rms=10;
			}else if(row==1){
				moved=10;
				better=8;
				worse=2;
				exact=110;
				rms=9;
			}else if(row==2){
				moved=8;
				better=7;
				worse=1;
				exact=108;
				rms=9.2;
			}else if(row==3){
				moved=20;
				better=15;
				worse=5;
				exact=115;
				rms=8.5;
			}else{
				final int delta=seed==58 ? -1 : seed==59 ? 7 : 0;
				if(cutoff==4){
					exact=120+delta;
					better=18;
					worse=2;
					rms=8;
				}
				if(cutoff==5){
					exact=113+delta;
					better=19;
					worse=0;
					rms=8.4;
				}
				if(cutoff==6){
					exact=112+delta;
					better=17;
					worse=2;
					rms=7;
				}
				if(cutoff==7){exact=120+delta;}
				if(outer && cutoff==9){
					exact=999;
					rms=1;
				}
			}
			values[5]=""+moved;
			values[6]=""+exact;
			values[9]=""+(100.0*exact/trials);
			values[11]="10";
			values[12]=""+rms;
			values[13]=""+(rms-10);
			values[14]=""+better;
			values[15]=""+worse;
			values[33]=""+(100.0*moved/trials);
			values[36]=""+(moved-better-worse);
			values[39]=""+(100.0*(better-worse)/trials);
			for(int i=0; i<40; i++){if(i>0){text.tab();}text.append(values[i]);}
			text.nl();
		}
		return text.toString();
	}
	/*--------------------------------------------------------------*/
	/*----------------          File Helpers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Invoke the actual production parser and table reader, not a copied algorithm. */
	private static void collect(final File root, final String name){
		new NcrnaGateCutoffSummary(new String[]{"in="+new File(root, "cells.tsv"), "outer="+new File(root, "outer"), "out="+new File(root, name)}).process();
	}
	/** Replace only a synthetic fixture receipt to exercise a negative case. */
	private static void write(final File file, final String text){
		try(Output out=new Output(new ByteStreamWriter(file.getPath(), true, false, false))){out.writer.print(text);}
	}
	/** Read the complete small output to check selection and literal values. */
	private static String read(final File file){
		final ByteBuilder result=new ByteBuilder();
		try(Input in=new Input(file.getPath())){for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){result.append(line).nl();}}
		return result.toString();
	}
	/** A rejected input is accepted only when the production tool throws loudly. */
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid saved-table evidence must fail");
	}
	private static void check(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
