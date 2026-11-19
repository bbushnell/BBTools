package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Exports exact original-locus exclusions for genome-level boundary validation.
 * Uses the ByteFile template's Parser/byte-reader/owned-writer lifecycle and the
 * existing strict population loader. Canonical and original FASTA paths join only
 * through the explicit corpus manifest; basenames, TaxIDs and sequence similarity
 * never stand in for assembly identity. No genome sequence payload or training
 * is read; original consensus names are loaded through the native library reader.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaLearningExclusions{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaLearningExclusions.class, false);
		try{
			new NcrnaLearningExclusions(pp.args).process();
			pp.outstream.println("NCRNA_LEARNING_EXCLUSIONS_COMPLETE exact_manifest_paths=true");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private NcrnaLearningExclusions(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("panel")){panel=b;}
			else if(a.equals("optionalpanel")){optionalPanel=b;}
			else if(a.equals("graph")){graph=b;}
			else if(a.equals("original")){original=b;}
			else if(a.equals("roster")){roster=b;}
			else if(a.equals("loci")){loci=b;}
			else{check(parser.parse(arg, a, b), "Unknown exclusion-export argument: "+arg);}
		}
		manifest=parser.in1;
		check(manifest!=null && parser.out1!=null && panel!=null
			&& graph!=null && original!=null && roster!=null && loci!=null,
			"Require in=corpus_manifest panel=required_FASTA_paths optionalpanel=other_FASTA_paths graph= original= roster= loci= out=fresh_directory");
		check(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Export complete source-qualified exclusions without sampling or overwrite");
		check(Tools.testInputFiles(false, true, manifest, graph, original, roster, loci)
			&& Tools.testInputFiles(false, true, panel, optionalPanel), "Unreadable source or panel binding");
		root=new File(parser.out1).getAbsoluteFile();
		check(!root.exists() && root.getParentFile()!=null && root.getParentFile().isDirectory(), "Exclusion output must be fresh beneath an existing parent");
	}

	private void process(){
		readManifest();
		final NcrnaLearningPopulation population=NcrnaLearningPopulation.loadForRefit(NcrnaSupportIO.read(graph),
			NcrnaSupport.loadLibrary(original).names, roster, loci, null);
		check(root.mkdir(), "Could not create exclusion output directory");
		final boolean[] heldout=new boolean[assemblyIds.size()];
		try(Output audit=new Output(path("panel_sources.tsv"))){
			audit.writer.println(LABEL);
			audit.writer.println("panel\trow1\tfasta\tassembly\tstatus");
			final ByteBuilder row=new ByteBuilder();
			for(final String panelFile:new String[]{panel, optionalPanel}){
				if(panelFile==null){continue;}
				int ordinal=0;
				try(Input input=new Input(panelFile)){
					for(byte[] line; (line=input.reader.nextLine())!=null;){
						final String fasta=new String(line, StandardCharsets.UTF_8);
						check(fasta.startsWith("/") && fasta.indexOf('\t')<0, "Panel inputs must be exact absolute FASTA paths");
						final int assembly=byPath.get(fasta);
						check(assembly>=0 || !panelFile.equals(panel), "A required-panel genome has no exact corpus alias: "+fasta);
						if(assembly>=0){heldout[assembly]=true;}
						else{unmatched++;}
						audit.writer.print(row.clear().append(panelFile).tab().append(++ordinal).tab().append(fasta).tab()
							.append(assembly<0 ? "." : assemblyIds.get(assembly)).tab()
							.append(assembly<0 ? "UNMATCHED_MANIFEST_PATH_ALIAS" : "EXCLUDED_ASSEMBLY").nl());
					}
				}
				check(ordinal>0, "Each evaluation panel must be nonempty");
			}
		}
		long excluded=0, eligibleBefore=0, eligibleAfter=0;
		try(Output output=new Output(path("exclude.txt"))){
			final ByteBuilder row=new ByteBuilder();
			for(int gene=0; gene<population.ids.length; gene++){
				final int assembly=byAssembly.get(population.assemblyName(gene));
				check(assembly>=0, "Original locus assembly is outside the bound working corpus: "+population.ids[gene]);
				if(population.eligible(gene)){eligibleBefore++;}
				if(heldout[assembly]){
					output.writer.print(row.clear().append(population.ids[gene]).nl());
					excluded++;
				}else if(population.eligible(gene)){eligibleAfter++;}
			}
		}
		long excludedAssemblies=0;
		for(final boolean selected:heldout){if(selected){excludedAssemblies++;}}
		check(excluded>0 && eligibleAfter>0 && eligibleAfter<=eligibleBefore, "Held-out comparison needs excluded loci and a nonempty retained teaching population");
		try(Output output=new Output(path("summary.tsv"))){
			output.writer.println(LABEL);
			output.writer.println("original_loci\texcluded_original_loci\tretained_original_loci\teligible_before\teligible_after\texcluded_manifest_assemblies\tpanel_paths_unmatched");
			output.writer.print(new ByteBuilder().append(population.ids.length).tab().append(excluded).tab()
				.append(population.ids.length-excluded).tab().append(eligibleBefore).tab().append(eligibleAfter)
				.tab().append(excludedAssemblies).tab().append(unmatched).nl());
		}
	}

	/** Only the manifest's two explicit path aliases identify an assembly. */
	private void readManifest(){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(manifest)){
			RrnaResourceIO.header(input.reader, MANIFEST_HEADER, manifest);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				p.set(line);
				check(p.terms()==12, "Working corpus manifest schema changed");
				final String id=p.parseString(0);
				check(!id.isEmpty() && !byAssembly.contains(id), "Corpus assembly IDs must be nonempty and unique");
				final int index=assemblyIds.size();
				byAssembly.put(id, index);
				assemblyIds.add(id);
				alias(p.parseString(3), index);
				alias(p.parseString(7), index);
			}
		}
		check(!assemblyIds.isEmpty(), "Working corpus manifest must contain assemblies");
	}
	private void alias(final String path, final int assembly){
		final int previous=byPath.get(path);
		check(path.startsWith("/") && (previous<0 || previous==assembly), "A manifest FASTA alias cannot identify different assemblies: "+path);
		byPath.put(path, assembly);
	}
	private String path(final String name){return new File(root, name).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	private String graph, original, roster, loci;
	private String panel, optionalPanel;
	private final String manifest;
	private final File root;
	private final ObjectIntMap<String> byPath=new ObjectIntMap<String>(String.class);
	private final ObjectIntMap<String> byAssembly=new ObjectIntMap<String>(String.class);
	private final ArrayList<String> assemblyIds=new ArrayList<String>();
	private long unmatched;
	private static final String LABEL="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result";
	private static final String MANIFEST_HEADER="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
}
