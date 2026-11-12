package prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Membership;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Selects ordered evaluation metadata without reading genome or annotation payloads.
 * Invoked through BuildNcrnaReference's native template-derived CLI. A declared
 * catalog source limits eligible inputs; the flat caller config selects and orders
 * them. Original FASTA aliases survive in the binder roster, canonical paths in
 * the grading manifest. No loci or masks are copied, inferred or reclassified.
 * The project recipe verifies registered identity before invoking this metadata
 * operation; a structurally valid manifest alone does not establish authorization.
 * One invocation owns all state and native readers/writers. The final reference
 * manifest is written only after metadata validation and other outputs close.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaEvaluationPanel{

	static void process(final NcrnaEvidenceCatalog catalog, final String config,
		final String source, final String reference, final String observed, final String prefix){
		assert(catalog!=null && catalog.size()>0) : "Panel selection requires the complete validated catalog";
		final NcrnaEvaluationPanel panel=new NcrnaEvaluationPanel(catalog, prefix);
		panel.select(config, source);
		panel.readReference(reference);
		panel.writeContigs(observed);
		panel.writeRoster();
		panel.writeSelection(config, source, reference);
		panel.writeReference();
	}

	private NcrnaEvaluationPanel(final NcrnaEvidenceCatalog catalog_, final String prefix_){
		catalog=catalog_;
		prefix=prefix_;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Selection           ----------------*/
	/*--------------------------------------------------------------*/

	/** A supplied path may alias one input, never two; repeated inputs fail even
	 * when their caller paths differ. Nested or comma-list configs are unsupported. */
	private void select(final String path, final String source){
		final ObjectMap<String, Selection> eligible=new ObjectMap<String, Selection>(String.class, Selection.class);
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			for(int j=0; j<genome.memberships(); j++){
				final Membership member=genome.membership(j);
				if(!member.source.equals(source)){continue;}
				final Selection selection=new Selection(genome, member);
				addPath(eligible, genome.fasta, selection);
				addPath(eligible, member.suppliedFasta, selection);
			}
		}
		require(eligible.size()>0, "Unknown or empty catalog source panel: "+source);
		try(ConfigInput config=new ConfigInput(path)){
			while(config.next()){
				final String key=config.key;
				if(key.equalsIgnoreCase("in") || key.equalsIgnoreCase("in1")){
					final String fasta=config.value;
					require(fasta.indexOf(',')<0, "Panel selection requires one literal FASTA per input record");
					final Selection selection=eligible.get(fasta);
					require(selection!=null, "Caller FASTA is not a member of sourcepanel: "+fasta);
					require(!selected.contains(selection.genome.id), "Repeated catalog input in caller configuration: "+selection.genome.id);
					selection.originalFasta=fasta;
					selected.put(selection.genome.id, selection);
					order.add(selection);
				}else{
					require(!key.equalsIgnoreCase("config") && !key.equalsIgnoreCase("in2"),
						"Panel preparation requires a flat unpaired caller input list");
				}
			}
		}
		require(!order.isEmpty(), "Caller configuration must declare inputs, including zero-call genomes");
	}

	private static void addPath(final ObjectMap<String, Selection> paths, final String path, final Selection selection){
		final Selection old=paths.get(path);
		require(old==null || old.genome==selection.genome
			&& old.member.sourceGenome.equals(selection.member.sourceGenome)
			&& old.member.assembly.equals(selection.member.assembly),
			"Ambiguous FASTA ownership or conflicting source metadata: "+path);
		if(old==null){paths.put(path, selection);}
	}

	/** Reads the complete manifest so a missing off-panel row cannot masquerade as
	 * the full registered population. Payload identity is the registration's job. */
	private void readReference(final String path){
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		final ObjectSet<String> files=new ObjectSet<String>(String.class);
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, GradeNcrnaReference.MANIFEST_HEADER, path);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==3, "Reference manifest requires input_id, fasta and gff");
				final String id=text(0), fasta=text(1), gff=text(2);
				final Genome genome=catalog.get(id);
				require(genome!=null && genome.fasta.equals(fasta) && ids.add(id),
					"Reference input must bind once to its exact catalog FASTA: "+id);
				require(new File(gff).isAbsolute() && files.add(gff), "Reference GFF paths must be absolute and unique");
				final Selection selection=selected.get(id);
				if(selection!=null){selection.gff=gff;}
			}
		}
		require(ids.size()==catalog.size(), "Reference manifest must cover the complete registered catalog");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Metadata output       ----------------*/
	/*--------------------------------------------------------------*/

	/** Streams saved contigs once. Only selected names are retained for the combined
	 * caller's global-uniqueness constraint; off-panel rows allocate no name strings. */
	private void writeContigs(final String path){
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		String previous=null;
		byte[] previousFasta=null;
		Selection selection=null;
		try(Input input=new Input(path); Output output=new Output(prefix+".contigs.tsv")){
			RrnaResourceIO.header(input.reader, WriteNcrnaReference.CONTIG_HEADER, path);
			output.writer.println(BindNcrnaCalls.PANEL_CONTIGS);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==4, "Saved contigs require four fields");
				if(previous==null || !parser.termEquals(previous, 0)){
					previous=text(0);
					selection=selected.get(previous);
					if(selection!=null){
						require(decoded(1).equals(selection.genome.fasta), "Observed FASTA differs from selected catalog input");
						parser.setBounds(1);
						previousFasta=Arrays.copyOfRange(row, parser.a(), parser.b());
					}
				}
				if(selection==null){continue;}
				require(sameField(1, previousFasta) || decoded(1).equals(selection.genome.fasta),
					"Observed FASTA differs from selected catalog input");
				final String name=decoded(2);
				final long length=parser.parseLong(3);
				require(!name.isEmpty() && names.add(name), "Repeated selected contig name would make combined calls ambiguous: "+name);
				require(length>0 && length<=Integer.MAX_VALUE, "Selected contig length must fit caller coordinates");
				selection.contigs++;
				selection.bases=Math.addExact(selection.bases, length);
				buffer.clear();
				raw(selection.genome.id).tab();
				raw(selection.member.assembly).tab();
				raw(name).tab().append(length).nl();
				output.writer.print(buffer);
			}
		}
		for(final Selection item:order){
			require(item.contigs==item.genome.records && item.bases==item.genome.bases,
				"Selected observed contigs/bases must cover the entire input: "+item.genome.id);
		}
	}

	/** Roster genome names are canonical input IDs, not guessed biological names.
	 * Unavailable domain/phylum/GFF/locus count remain NA rather than zero truth. */
	private void writeRoster(){
		try(Output output=new Output(prefix+".genomes.tsv")){
			output.writer.println(NcrnaCalibrationPrep.ROSTER_HEADER);
			int index=0;
			for(final Selection item:order){
				buffer.clear().append(++index).tab();
				raw(item.genome.id).tab();
				raw(item.member.assembly).tab();
				raw(item.member.taxid).append("\tNA\tNA\tNA\t");
				raw(item.originalFasta).append("\tNA\t").append(item.bases).tab().append(item.contigs).nl();
				output.writer.print(buffer);
			}
		}
	}

	/** Records selection, not a claim that the Java metadata reader verified hashes. */
	private void writeSelection(final String config, final String source, final String reference){
		try(Output output=new Output(prefix+".selection.tsv")){
			output.writer.println("input_id\tsource_panel\tcaller_fasta\tcanonical_fasta\treference_gff\tcaller_config\treference_manifest");
			for(final Selection item:order){
				assert(item.gff!=null) : "Complete reference-catalog validation precedes selection output";
				buffer.clear();
				raw(item.genome.id).tab();
				raw(source).tab();
				raw(item.originalFasta).tab();
				raw(item.genome.fasta).tab();
				raw(item.gff).tab();
				raw(config).tab();
				raw(reference).nl();
				output.writer.print(buffer);
			}
		}
	}

	private void writeReference(){
		assert(!order.isEmpty()) : "A completed panel manifest must contain every selected input";
		try(Output output=new Output(prefix+".reference.tsv")){
			output.writer.println(GradeNcrnaReference.MANIFEST_HEADER);
			for(final Selection item:order){
				buffer.clear();
				raw(item.genome.id).tab();
				raw(item.genome.fasta).tab();
				raw(item.gff).nl();
				output.writer.print(buffer);
			}
		}
	}

	/** Literal TSV framing matches the binder/grade manifest contract; UTF-8 is
	 * preserved rather than truncated by ByteBuilder's ASCII String overload. */
	private ByteBuilder raw(final String value){
		require(value!=null && value.indexOf('\t')<0 && value.indexOf('\r')<0 && value.indexOf('\n')<0,
			"Literal panel metadata cannot contain TSV record separators");
		for(int i=0; i<value.length(); i++){
			if(value.charAt(i)>127){return buffer.append(value.getBytes(StandardCharsets.UTF_8));}
		}
		return buffer.append(value);
	}
	private String text(final int field){
		parser.setBounds(field);
		final String result=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(!result.isEmpty() && result.indexOf('\ufffd')<0, "Panel metadata needs nonempty lossless UTF-8 fields");
		return result;
	}
	private String decoded(final int field){
		parser.setBounds(field);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	/** Compares a repeated byte field without allocating one path String per contig. */
	private boolean sameField(final int field, final byte[] expected){
		assert(expected!=null) : "The first selected contig binds its input path before repeated rows";
		parser.setBounds(field);
		if(parser.b()-parser.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){
			if(parser.line()[parser.a()+i]!=expected[i]){return false;}
		}
		return true;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}

	private static final class Selection{
		Selection(final Genome genome_, final Membership member_){
			genome=genome_;
			member=member_;
		}
		final Genome genome;
		final Membership member;
		String originalFasta, gff;
		long contigs, bases;
	}
	private final NcrnaEvidenceCatalog catalog;
	private final String prefix;
	private final ArrayList<Selection> order=new ArrayList<Selection>();
	private final ObjectMap<String, Selection> selected=new ObjectMap<String, Selection>(String.class, Selection.class);
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder buffer=new ByteBuilder();
	static final String[] SUFFIXES={".genomes.tsv", ".contigs.tsv", ".selection.tsv", ".reference.tsv"};
}
