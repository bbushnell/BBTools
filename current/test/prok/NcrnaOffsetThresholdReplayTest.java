package test.prok;

import java.io.File;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.EvaluateNcrnaRealEndpoints;
import prok.NcrnaEndpointFoldInputs;
import prok.NcrnaOffsetThresholdReplay;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Literal signed offsets, physical guards and an outer population deliberately
 * preferring the wrong monitoring choice detect threshold/selection leakage.
 * No real genomes, models, alignments or fits.
 * @author Brian, Ganyu
 */
public final class NcrnaOffsetThresholdReplayTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaOffsetThresholdReplayTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown test option");
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(),
				"Supply out=fresh_fixture_directory");
			policies();
			fixture(new File(parser.out1));
			pp.outstream.println("NCRNA_OFFSET_THRESHOLD_REPLAY_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void policies(){
		final int[] ends=new int[2];
		check(NcrnaOffsetThresholdReplay.apply(30, 149, 250, 0, 109, 200, 5, -3, 3, 2, ends)==0
			&& ends[0]==25 && ends[1]==149, "Strict threshold excludes magnitude equal to T and subtracts signed medians");
		check(NcrnaOffsetThresholdReplay.apply(30, 149, 250, 0, 109, 200, 5, -3, 2, 2, ends)==0
			&& ends[0]==25 && ends[1]==152, "Negative stop median extends in gene sense");
		check(NcrnaOffsetThresholdReplay.apply(0, 149, 250, 1, 109, 200, 5, 3, 0, 2, ends)==0
			&& ends[0]==0 && ends[1]==146, "Physical fixed ends cannot move even when supported");
		check(NcrnaOffsetThresholdReplay.apply(30, 149, 250, 0, 109, 200, -10, 10, 0, 2, ends)==3
			&& ends[0]==30 && ends[1]==149, "Joint shortening below the family floor restores both raw ends");
		check(NcrnaOffsetThresholdReplay.apply(30, 149, 250, 0, 109, 200, -10, 10, 0, 0, ends)==0
			&& ends[0]==40 && ends[1]==149, "The same one-end move can be valid alone with the other end raw");
		check(NcrnaOffsetThresholdReplay.apply(2, 149, 250, 0, 109, 200, 5, 3, 0, 2, ends)==2
			&& ends[0]==2 && ends[1]==149, "An outside-contig center holds rather than clipping");
		check(NcrnaOffsetThresholdReplay.apply(30, 149, 250, 0, 109, 200, 5, 3, Integer.MAX_VALUE, 2, ends)==0
			&& ends[0]==30 && ends[1]==149, "Never-shift control preserves both raw coordinates");
	}

	private static void fixture(final File root){
		final String model="cobalamin_RF00174_fixture";
		final String reference=path(root, "reference.fa");
		write(reference, new ByteBuilder().append('>').append(model).append("\nACGTACGTACGT\n"));
		final ByteBuilder folds=new ByteBuilder(NcrnaEndpointFoldInputs.ASSEMBLY_HEADER).nl();
		final ByteBuilder train=new ByteBuilder(EvaluateNcrnaRealEndpoints.GEOMETRY_HEADER).nl();
		for(int fold=0; fold<5; fold++){
			folds.append('a').append(fold).append("\t/fixture/").append(fold).append(".fa\tINTERNAL_FOLD\t").append(fold).nl();
			for(int i=0; i<3; i++){
				// Held fold4 is already exact; the other folds support a nonzero
				// median. It must reject movement even when the outer set favors it.
				geometry(train, fold*3+i+1, "a"+fold, "ref"+i, model, fold==4 ? 0 : 5, fold==4 ? 0 : 3);
			}
		}
		folds.append("outer\t/fixture/outer.fa\tOUTER_TEST\t-1\n");
		train.append("#COMPLETE\t15\n");
		final ByteBuilder outer=new ByteBuilder(EvaluateNcrnaRealEndpoints.GEOMETRY_HEADER).nl();
		geometry(outer, 1, "outer", "r", model, 5, 3);
		outer.append("#COMPLETE\t1\n");
		write(path(root, "folds.tsv"), folds);
		write(path(root, "train.tsv"), train);
		write(path(root, "new.tsv"), outer);
		write(path(root, "old.tsv"), outer);
		final ByteBuilder manifest=new ByteBuilder(NcrnaOffsetThresholdReplay.DATASET_HEADER).nl();
		manifest.append("train\tTRAIN\t").append(path(root, "train.tsv")).append("\t15\n")
			.append("new\tOUTER\t").append(path(root, "new.tsv")).append("\t1\n")
			.append("old\tOUTER\t").append(path(root, "old.tsv")).append("\t1\n");
		write(path(root, "datasets.tsv"), manifest);
		final String[] args={"in="+path(root, "datasets.tsv"), "folds="+path(root, "folds.tsv"),
			"ref="+reference, "member=RF00174", "out="+path(root, "result")};
		new NcrnaOffsetThresholdReplay(args).process();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path(root, "result/fold4.offsets.tsv"))){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.parseInt(2)==12 && p.parseInt(3)==10 && p.parseInt(4)==5 && p.parseInt(5)==3,
				"Fold4 offsets must exclude its three exact rows and retain the twelve other real calls");
		}
		int selected=0;
		try(Input in=new Input(path(root, "result/thresholds.tsv"))){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				final boolean chosen=p.termEquals("1", 4);
				check(chosen==(p.termEquals("5", 1) || p.termEquals("never", 1)),
					"Only zero-movement tied monitor optima may be selected; outer improvement cannot select threshold0");
				if(chosen){selected++;}
			}
		}
		check(selected==2, "Retain both tied monitor optima instead of using outer truth to break ties");
		int selectedOuter=0;
		try(Input in=new Input(path(root, "result/outer.tsv"))){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==36 && p.parseInt(6)==1, "Every outer policy retains the single physical call");
				if(p.termEquals("1", 3)){
					check(p.parseInt(19)==0 && p.parseInt(7)==0, "Frozen selected thresholds hold even though outer truth would reward a move");
					selectedOuter++;
				}
			}
		}
		check(selectedOuter==16, "Two selected thresholds times two panels times four single/joint end rows");
		int audited=0;
		try(Input in=new Input(path(root, "result/outer_decisions.tsv"))){
			in.reader.nextLine();
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==23 && p.termEquals("input_outer", 5) && p.termEquals("outer", 6)
					&& p.termEquals("r", 7) && p.parseInt(13)==15 && p.parseInt(14)==5 && p.parseInt(15)==3,
					"Outer audit preserves exact identities and full-training support/medians");
				final boolean never=!p.termEquals("0", 2);
				check(p.parseInt(21)==(never || p.termEquals("SINGLE_3PRIME", 4) ? 50 : 45)
					&& p.parseInt(22)==(never || p.termEquals("SINGLE_5PRIME", 4) ? 179 : 176),
					"Audit coordinates retain the literal separate-end and paired policies");
				audited++;
			}
		}
		check(audited==18, "Every selected/control threshold has three scopes on both outer panels");
		residuals(root, args, manifest);
		// A support assembly absent from the mapping must never default into
		// training. Change just that identity and demand failure before output.
		write(path(root, "bad_train.tsv"), new ByteBuilder(train.toString().replace("\ta0\t", "\tunknown\t")));
		write(path(root, "bad_datasets.tsv"), new ByteBuilder(manifest.toString().replace(path(root, "train.tsv"), path(root, "bad_train.tsv"))));
		args[0]="in="+path(root, "bad_datasets.tsv");
		args[4]="out="+path(root, "bad_unknown");
		boolean failed=false;
		try{new NcrnaOffsetThresholdReplay(args).process();}
		catch(IllegalArgumentException expected){failed=expected.getMessage().contains("Unknown canonical assembly");}
		check(failed && !new File(path(root, "bad_unknown")).exists(), "Unknown support assemblies fail before results are created");
	}

	/** Residual export must use accepted medians, exclude the held fold and
	 * avoid opening outer truth; deliberately nonexistent OUTER files prove it. */
	private static void residuals(final File root, final String[] original, final ByteBuilder manifest){
		final String datasets=path(root, "residual_datasets.tsv");
		write(datasets, new ByteBuilder(manifest.toString()
			.replace(path(root, "new.tsv"), path(root, "outer_new_must_not_be_opened.tsv"))
			.replace(path(root, "old.tsv"), path(root, "outer_old_must_not_be_opened.tsv"))));
		final String[] args=Arrays.copyOf(original, original.length+2);
		args[0]="in="+datasets;
		args[4]="out="+path(root, "residuals");
		args[original.length]="mode=residuals";
		args[original.length+1]="offsets="+path(root, "result");
		new NcrnaOffsetThresholdReplay(args).process();
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(path(root, "residuals/fold4.histogram.tsv"))){
			in.reader.nextLine();
			for(final String end:new String[]{"5prime", "3prime"}){
				p.set(in.reader.nextLine());
				check(p.terms()==5 && p.termEquals("RF00174", 0) && p.termEquals(end, 1)
					&& p.termEquals("signed", 2) && p.parseInt(3)==0 && p.parseInt(4)==12,
					"The twelve non-held witnesses become zero residuals after their own model's accepted offset");
			}
			check(in.reader.nextLine()==null, "Held-fold truths cannot contribute residual weight");
		}
		try(Input in=new Input(path(root, "residuals/fold0.histogram.tsv"))){
			in.reader.nextLine();
			for(final int offset:new int[]{-5, 0, -3, 0}){
				p.set(in.reader.nextLine());
				check(p.parseInt(3)==offset && p.parseInt(4)==(offset==0 ? 9 : 3),
					"Residual counts preserve signed overcorrections rather than censoring hard training examples");
			}
			check(in.reader.nextLine()==null, "Each retained call contributes exactly one count per end");
		}
		check(!new File(path(root, "residuals/thresholds.tsv")).exists()
			&& !new File(path(root, "residuals/outer.tsv")).exists(), "Residual export never repeats selection or evaluation");
	}

	private static void geometry(final ByteBuilder row, final int ordinal, final String assembly,
			final String ref, final String model, final int error5, final int error3){
		row.append(ordinal).append("\tinput_").append(assembly).tab().append(assembly).tab().append(ref).append("\tRF00174\t")
			.append(model).tab().append(ref).append("\t+\t300\t0\t50\t179\t")
			.append(50-error5).tab().append(179-error3).append("\t109\t\n");
	}
	private static String path(final File root, final String name){return new File(root, name).getAbsolutePath();}
	private static void write(final String path, final ByteBuilder data){try(Output out=new Output(path)){out.writer.print(data);}}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
