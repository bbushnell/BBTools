package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SelectNcrnaBenchmarkInputs;
import shared.Shared;

/** Finite metadata-only selection oracles; opens no genome payloads.
 * Checks proportional quotas, explicit unknowns, deterministic taxonomy reorder,
 * complete-population identity, caller declarations and malformed populations.
 * @author Brian Bushnell, Raiden
 */
public final class SelectNcrnaBenchmarkInputsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, SelectNcrnaBenchmarkInputsTest.class, false);
		final Parser parser=new Parser();
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
			if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		require(parser.out1!=null, "Supply out=fresh_fixture_directory");
		final File root=new File(parser.out1);
		require(!root.exists(), "Fixture directory must be fresh");
		Files.createDirectory(root.toPath());
		try{test(root);}
		finally{Shared.closeStream(pp.outstream);}
		System.err.println("NCRNA_BENCHMARK_INPUTS_TEST_PASS");
	}

	private static void test(final File root) throws Exception{
		final String census=path(root, "genomes.tsv"), memberships=path(root, "memberships.tsv"), taxonomy=path(root, "taxonomy.tsv");
		final ArrayList<String> tax=new ArrayList<String>();
		tax.add("input_id\tfasta\tdomain\tphylum\tcorpus_inputs\tcategory\ttaxonomy_status");
		try(Output genomes=new Output(census); Output members=new Output(memberships)){
			text(genomes, NcrnaEvidenceCatalog.CENSUS_HEADER);
			text(members, NcrnaEvidenceCatalog.MEMBERSHIP_HEADER);
			for(int i=1; i<=14; i++){
				final String fasta=fasta(i), group=i<=8 ? "A" : i<=12 ? "B" : "UNKNOWN";
				text(genomes, i+"\t"+fasta+"\t100\t"+(1000L*i)+"\t"+i);
				text(members, "control\t"+i+"\t"+i+"\t"+fasta+"\tgenome"+i+"\tassembly"+i+"\tNA\tprokaryotic\tNA\tNA");
				tax.add(i+"\t"+fasta+"\t"+(i<=12 ? "bacteria" : "UNKNOWN")+"\t"+group+"\t14\tUNUSED\tEXPLICIT");
			}
		}
		write(taxonomy, tax);
		run(census, memberships, taxonomy, path(root, "first"), "7,14");
		final String selected=path(root, "first.7.selection.tsv");
		final int[] counts=new int[3];
		long selectedRecords=0;
		final LineParser1 parser=new LineParser1('\t');
		try(Input in=new Input(selected)){
			require(in.reader.nextLine()!=null, "Selection header must be present");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==8, "Selection must retain complete workload and taxonomy metadata");
				final int id=parser.parseInt(0);
				counts[id<=8 ? 0 : id<=12 ? 1 : 2]++;
				selectedRecords+=id;
			}
		}
		require(Arrays.equals(counts, new int[]{4, 2, 1}),
			"Largest-remainder allocation of 7 slots among populations 8/4/2 must be 4/2/1, retaining unknowns");
		final ArrayList<String> subsetTotals=lines(path(root, "first.7.totals.tsv"));
		require(subsetTotals.size()==2 && subsetTotals.get(1).equals(
			"14\t7\t105000\t"+(1000*selectedRecords)+"\t105\t"+selectedRecords),
			"Subset totals must follow the selected fixture IDs: each has id records and 1000 times id bases");
		final ArrayList<String> all=lines(path(root, "first.14.inputs.txt"));
		require(all.size()==14, "The full-size selection must retain every original input");
		for(int i=0; i<all.size(); i++){require(all.get(i).equals(fasta(i+1)), "Full-size selection must preserve original input order and UTF-8");}
		final ArrayList<String> flat=lines(path(root, "first.14.flat.config"));
		for(int i=0; i<all.size(); i++){require(flat.get(i).equals("in="+all.get(i)), "Flat declarations must preserve exact paths");}
		require(lines(path(root, "first.14.inputs.config")).get(0).equals("in="+String.join(",", all)),
			"One config record must preserve the complete comma-separated input population");
		final long[] totals=new long[6];
		try(Input in=new Input(path(root, "first.14.composition.tsv"))){
			require(in.reader.nextLine()!=null, "Composition needs its declared schema");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				for(int i=0; i<totals.length; i++){totals[i]+=parser.parseLong(i+4);}
			}
		}
		require(Arrays.equals(totals, new long[]{14, 14, 105000, 105000, 105, 105}),
			"Joint-bin accounting must conserve the literal full census dimensions, not merely output row counts");
		final ArrayList<String> aggregate=lines(path(root, "first.14.totals.tsv"));
		require(aggregate.size()==2 && aggregate.get(0).equals(
			"corpus_inputs\tselected_inputs\tcorpus_bases\tselected_bases\tcorpus_records\tselected_records")
			&& aggregate.get(1).equals("14\t14\t105000\t105000\t105\t105"),
			"Aggregate totals must report the independent literal full census dimensions in the declared units");

		final ArrayList<String> reversed=new ArrayList<String>();
		reversed.add(tax.get(0));
		for(int i=tax.size()-1; i>0; i--){reversed.add(tax.get(i));}
		final String reordered=path(root, "reordered.tsv");
		write(reordered, reversed);
		run(census, memberships, reordered, path(root, "second"), "7,14");
		for(final int size:new int[]{7, 14}){
			for(final String suffix:new String[]{"inputs.txt", "inputs.config", "flat.config", "selection.tsv", "composition.tsv", "totals.tsv"}){
				require(Arrays.equals(Files.readAllBytes(new File(path(root, "first."+size+"."+suffix)).toPath()),
					Files.readAllBytes(new File(path(root, "second."+size+"."+suffix)).toPath())),
					"Reordering taxonomy must not change selection, workload accounting or original input order");
			}
		}
		final ArrayList<String> duplicate=new ArrayList<String>(tax);
		duplicate.add(tax.get(1));
		write(path(root, "duplicate.tsv"), duplicate);
		reject(census, memberships, path(root, "duplicate.tsv"), path(root, "duplicate"), "7", "once");
		final ArrayList<String> missing=new ArrayList<String>(tax);
		missing.remove(missing.size()-1);
		write(path(root, "missing.tsv"), missing);
		reject(census, memberships, path(root, "missing.tsv"), path(root, "missing"), "7", "Every census input");
		reject(census, memberships, taxonomy, path(root, "oversized"), "15", "cannot exceed");
		reject(census, memberships, taxonomy, path(root, "toosmall"), "2", "every domain/phylum");
	}

	private static void run(final String census, final String memberships, final String taxonomy, final String prefix, final String sizes){
		SelectNcrnaBenchmarkInputs.main(new String[]{"in="+census, "memberships="+memberships,
			"taxonomy="+taxonomy, "sizes="+sizes, "out="+prefix, "seed=7", "t=1"});
	}
	private static void reject(final String census, final String memberships, final String taxonomy, final String prefix, final String sizes, final String diagnostic){
		try{run(census, memberships, taxonomy, prefix, sizes);}
		catch(IllegalArgumentException failure){
			require(failure.getMessage().contains(diagnostic), "Wrong invalid-population diagnostic: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Invalid selection was accepted: "+diagnostic);
	}
	private static String path(final File root, final String name){return new File(root, name).toString();}
	private static String fasta(final int id){return "/genomes/α="+id+".fa";}
	private static void text(final Output out, final String line){out.writer.println(line.getBytes(StandardCharsets.UTF_8));}
	private static ArrayList<String> lines(final String path){
		final ArrayList<String> result=new ArrayList<String>();
		try(Input in=new Input(path)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){result.add(new String(row, StandardCharsets.UTF_8));}
		}
		return result;
	}
	private static void write(final String path, final ArrayList<String> rows){
		try(Output out=new Output(path)){for(final String row:rows){text(out, row);}}
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new AssertionError(message);}}
}
