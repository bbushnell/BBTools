package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEndpointFoldInputs;
import prok.NcrnaTextIO.Input;
import shared.Shared;
import structures.ByteBuilder;

/** Five-fold assembly conservation, deterministic assignment and outer leakage controls.
 * Uses the canonical original-population fixture; no genome or model computation.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaEndpointFoldInputsTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEndpointFoldInputsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, pair[0].toLowerCase(Locale.ROOT), pair.length>1 ? pair[1] : null), "Unknown fixture option");
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.maxReads<0
				&& parser.samplerate==1 && !parser.overwrite && !parser.append && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			run(new File(parser.out1));
			pp.outstream.println("NCRNA_ENDPOINT_FOLD_INPUTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String manifest=path(root, "manifest.tsv");
		final String roster=path(root, "roster.tsv");
		final String loci=path(root, "loci.tsv");
		final String graph=path(root, "graph.tsv");
		final String original=path(root, "original.fa");
		final String exclude=path(root, "exclude.txt");
		final String audit=path(root, "audit.tsv");
		write(roster, NcrnaLearningPopulationTest.roster());
		write(loci, NcrnaLearningPopulationTest.loci());
		write(graph, NcrnaLearningPopulationTest.graph());
		write(original, ">common\nAAAA\n>rare\nCCCC\n");
		write(exclude, "gene2\ngene3\n");
		final ByteBuilder corpus=new ByteBuilder(NcrnaEndpointFoldInputs.MANIFEST_HEADER).nl();
		for(int i=0; i<=20; i++){
			corpus.append("assembly").append(i).tab().append("genome").append(i)
				.append("\tfixture\t/canonical/").append(i).append(".fa\t/gff\t.\t.\t/alias/")
				.append(i).append(".fa\t/source.gff\tfasta\tncrna\tretained\n");
		}
		write(manifest, corpus.toString());
		final String auditText=NcrnaEndpointFoldInputs.AUDIT_LABEL+'\n'+NcrnaEndpointFoldInputs.AUDIT_HEADER
			+"\n/outer1000\t1\t/canonical/2.fa\tassembly2\tEXCLUDED_ASSEMBLY\n"
			+"/outer300\t1\t/alias/3.fa\tassembly3\tEXCLUDED_ASSEMBLY\n";
		write(audit, auditText);
		final String[] args={"in="+manifest, "out="+path(root, "first"), "original="+original,
			"graph="+graph, "roster="+roster, "loci="+loci, "exclude="+exclude, "audit="+audit, "seed=20261117"};
		new NcrnaEndpointFoldInputs(args).process();
		args[1]="out="+path(root, "second");
		new NcrnaEndpointFoldInputs(args).process();
		final int[] assigned=new int[21];
		Arrays.fill(assigned, -2);
		try(Input first=new Input(path(new File(root, "first"), "loci.tsv"));
				Input second=new Input(path(new File(root, "second"), "loci.tsv"))){
			check(Arrays.equals(first.reader.nextLine(), second.reader.nextLine()), "Deterministic mapping headers");
			final LineParser1 p=new LineParser1('\t');
			int count=0;
			for(byte[] line=first.reader.nextLine(); line!=null; line=first.reader.nextLine()){
				check(Arrays.equals(line, second.reader.nextLine()), "Same seed must reproduce exact assembly assignment independent of output path");
				p.set(line);
				check(p.terms()==6 && p.termEquals("gene"+count, 0), "Every original stays on the stable locus axis");
				assigned[count]=p.parseInt(4);
				check((count==2 || count==3) ? p.termEquals("OUTER_TEST", 3) && assigned[count]==-1
					: p.termEquals("INTERNAL_FOLD", 3) && assigned[count]>=0 && assigned[count]<5,
					"Outer genomes never acquire an internal fold");
				count++;
			}
			check(count==21 && second.reader.nextLine()==null && assigned[0]==assigned[20],
				"All originals and the partial sibling retain indivisible assembly membership");
		}
		final int[] heldCounts=new int[21];
		for(int fold=0; fold<5; fold++){
			final File folder=new File(root, "first/fold"+fold);
			final boolean[] seen=new boolean[21];
			try(Input input=new Input(path(folder, "exclude.txt"))){
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					final String id=new String(line, StandardCharsets.US_ASCII);
					final int gene=Integer.parseInt(id.substring(4));
					check(!seen[gene] && (assigned[gene]==-1 || assigned[gene]==fold), "Only outer and this fold's exact loci may be excluded");
					seen[gene]=true;
					heldCounts[gene]++;
				}
			}
			for(int gene=0; gene<21; gene++){
				check(seen[gene]==(assigned[gene]==-1 || assigned[gene]==fold), "No fold-held or outer locus may leak into table/NN teaching");
			}
			try(Input input=new Input(path(folder, "panel_sources.tsv"))){
				check(new String(input.reader.nextLine(), StandardCharsets.US_ASCII).equals(NcrnaEndpointFoldInputs.AUDIT_LABEL), "Audit label is native-consumable");
				input.reader.nextLine();
				check(new String(input.reader.nextLine(), StandardCharsets.US_ASCII).equals("/outer1000\t1\t/canonical/2.fa\tassembly2\tEXCLUDED_ASSEMBLY")
					&& new String(input.reader.nextLine(), StandardCharsets.US_ASCII).equals("/outer300\t1\t/alias/3.fa\tassembly3\tEXCLUDED_ASSEMBLY"),
					"Both outer audits remain intact before fold exclusions");
				final LineParser1 p=new LineParser1('\t');
				int ordinal=0;
				for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
					p.set(line);
					final int assembly=Integer.parseInt(p.parseString(3).substring(8));
					check(p.terms()==5 && p.parseInt(1)==++ordinal && assigned[assembly]==fold,
						"Histogram exclusions use the same fold and dense panel ordinals");
				}
				check(ordinal==3 || ordinal==4, "Eighteen internal assemblies divide into five nearly equal folds");
			}
		}
		for(int gene=0; gene<21; gene++){check(heldCounts[gene]==(assigned[gene]==-1 ? 5 : 1), "Outer loci excluded five times; each internal locus held once");}
		split(root, assigned);
		final String[] allArgs=Arrays.copyOf(args, args.length+1);
		allArgs[1]="out="+path(root, "all_assemblies");
		allArgs[args.length]="allassemblies=t";
		new NcrnaEndpointFoldInputs(allArgs).process();
		final int[] allFolds=new int[21];
		Arrays.fill(allFolds, -2);
		try(Input input=new Input(path(new File(root, "all_assemblies"), "assemblies.tsv"))){
			check(new String(input.reader.nextLine(), StandardCharsets.US_ASCII).equals(NcrnaEndpointFoldInputs.ASSEMBLY_HEADER),
				"The complete assembly map has a declared native schema");
			final LineParser1 p=new LineParser1('\t');
			int count=0;
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				final int assembly=Integer.parseInt(p.parseString(0).substring(8));
				check(p.terms()==4 && allFolds[assembly]==-2, "Each corpus assembly is assigned exactly once");
				allFolds[assembly]=p.parseInt(3);
				check(assembly==2 || assembly==3 ? allFolds[assembly]==-1 && p.termEquals("OUTER_TEST", 2)
					: allFolds[assembly]>=0 && allFolds[assembly]<5 && p.termEquals("INTERNAL_FOLD", 2),
					"Outer exclusion includes exactly the bound assemblies, even when unused for teaching");
				count++;
			}
			check(count==21 && allFolds[20]>=0, "The assembly without any teaching locus still needs a support fold");
		}
		try(Input input=new Input(path(new File(root, "all_assemblies"), "loci.tsv"))){
			input.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				check(p.parseInt(4)==allFolds[Integer.parseInt(p.parseString(2).substring(8))],
					"Teaching-locus exclusions and support witnesses must consume the same whole-assembly fold");
			}
		}
		write(audit, auditText.replace("/alias/3.fa", "/alias/4.fa"));
		args[1]="out="+path(root, "bad_alias");
		fails(()->new NcrnaEndpointFoldInputs(args).process());
		write(audit, auditText);
		write(exclude, "gene2\n");
		args[1]="out="+path(root, "bad_exclusion");
		fails(()->new NcrnaEndpointFoldInputs(args).process());
	}

	/** One outer plus one held row must become exactly one explicitly FOLD_HELD row. */
	private static void split(final File root, final int[] assigned){
		int held=-1;
		int train=-1;
		for(int gene=0; gene<20; gene++){
			if(assigned[gene]==0){held=gene;}
			if(assigned[gene]>0){train=gene;}
		}
		check(held>=0 && train>=0, "Literal split needs held/training assemblies");
		final File input=new File(root, "vectors");
		check(input.mkdir(), "Fresh fold vector fixture");
		for(final String set:new String[]{"train", "monitor", "evaluation"}){
			for(final String end:new String[]{"5prime", "3prime"}){
				final String stem="r10."+set+'.'+end;
				final boolean evaluation=set.equals("evaluation");
				write(path(input, stem+".tsv"), "#dims\t69\t22\t0\n"+(evaluation ? vector(2)+vector(held)+vector(2)+vector(held) : vector(train)));
				write(path(input, stem+".rows.tsv"), prok.NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'
					+(evaluation ? row(2, set)+row(held, set)+row(2, set, 1)+row(held, set, 1) : row(train, set)));
			}
		}
		final String[] args={"mode=split", "in="+input, "out="+path(root, "split_good"),
			"mapping="+path(new File(root, "first"), "loci.tsv"), "fold=0", "samples=2", "expected=1,1,1,1,4,4"};
		new NcrnaEndpointFoldInputs(args).process();
		for(final String end:new String[]{"5prime", "3prime"}){
			try(Input matrix=new Input(path(new File(root, "split_good"), "r10.fold_held."+end+".tsv.gz"));
					Input rows=new Input(path(new File(root, "split_good"), "r10.fold_held."+end+".rows.tsv"))){
				matrix.reader.nextLine();
				rows.reader.nextLine();
				check(new String(matrix.reader.nextLine(), StandardCharsets.US_ASCII).equals(vector(held).trim())
					&& new String(rows.reader.nextLine(), StandardCharsets.US_ASCII).equals(row(held, "FOLD_HELD").trim())
					&& new String(matrix.reader.nextLine(), StandardCharsets.US_ASCII).equals(vector(held).trim())
					&& new String(rows.reader.nextLine(), StandardCharsets.US_ASCII).equals(row(held, "FOLD_HELD", 1).trim())
					&& matrix.reader.nextLine()==null && rows.reader.nextLine()==null,
					"Outer vectors are discarded; held vector bytes survive and role is explicit");
			}
		}
		try(Input summary=new Input(path(new File(root, "split_good"), "split.tsv"))){
			summary.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			int records=0;
			for(byte[] line=summary.reader.nextLine(); line!=null; line=summary.reader.nextLine()){
				p.set(line);
				final boolean evaluation=p.termEquals("evaluation", 1);
				check(p.terms()==8 && p.parseLong(3)==(evaluation ? 4 : 1) && p.parseLong(4)==(evaluation ? 0 : 1)
					&& p.parseLong(5)==(evaluation ? 2 : 0) && p.parseLong(6)==(evaluation ? 2 : 0)
					&& p.parseLong(7)==(evaluation ? 2 : 0), "Every input row has one disjoint fold disposition");
				records++;
			}
			check(records==6, "Validate both ends in all three input splits");
		}
		// Both files agree on two valid rows and retain both cohorts, but the
		// bound producer emitted four. Paired EOF must not hide this truncation.
		final String evaluationMatrix=path(input, "r10.evaluation.5prime.tsv");
		final String evaluationRows=path(input, "r10.evaluation.5prime.rows.tsv");
		write(evaluationMatrix, "#dims\t69\t22\t0\n"+vector(2)+vector(held));
		write(evaluationRows, prok.NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(2, "evaluation")+row(held, "evaluation"));
		args[2]="out="+path(root, "paired_truncation");
		boolean rejected=false;
		try{new NcrnaEndpointFoldInputs(args).process();}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("expected=4 actual=2");}
		check(rejected, "Equal paired truncation must fail against producer counts, even with nonempty outer and held cohorts");
		final String trainingRows=path(input, "r10.train.5prime.rows.tsv");
		write(trainingRows, prok.NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(2, "train"));
		args[2]="out="+path(root, "outer_leak");
		fails(()->new NcrnaEndpointFoldInputs(args).process());
		write(trainingRows, prok.NcrnaPositionalFeatureVariants.ROW_HEADER+'\n'+row(held, "train"));
		args[2]="out="+path(root, "held_leak");
		fails(()->new NcrnaEndpointFoldInputs(args).process());
	}
	private static String vector(final int gene){
		final ByteBuilder row=new ByteBuilder().append(gene).tab();
		for(int i=1; i<69; i++){row.append("0\t");}
		for(int i=0; i<22; i++){
			if(i>0){row.tab();}
			row.append(i==10 ? 1 : 0);
		}
		return row.nl().toString();
	}
	private static String row(final int gene, final String split){
		return row(gene, split, 0);
	}
	private static String row(final int gene, final String split, final int sample){
		return "gene"+gene+"\tgenome"+gene+"\tassembly"+gene+'\t'+split+"\t0\tcommon\t"+sample+"\t20\t119\t20\t119\n";
	}
	private static String path(final File root, final String file){return new File(root, file).getPath();}
	private static void write(final String path, final String text){NcrnaLearningPopulationTest.write(path, text);}
	private static void fails(final Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Inconsistent fold binding must fail loudly");
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
}
