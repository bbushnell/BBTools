package prok;

import java.util.ArrayList;

import map.ObjectMap;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaEvidenceCatalog.Membership;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;

/** Counts saved panel memberships separately from distinct input IDs.
 * Does not open genomes or infer organism equivalence from assembly or TaxID.
 * Source names are the catalog's declared panel labels, not biological families;
 * assigning a shared panel to particular callers remains explicit project metadata.
 * One thread owns accumulation; the supplied catalog is borrowed read-only.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferencePanels{

	static void process(final NcrnaEvidenceCatalog catalog, final String prefix){
		assert(catalog!=null && catalog.size()>0) : "Panel accounting requires the validated input/membership catalog";
		final ObjectMap<String,Counts> sources=new ObjectMap<String,Counts>(String.class, Counts.class);
		final ArrayList<Counts> order=new ArrayList<Counts>();
		long memberships=0;
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			for(int j=0; j<genome.memberships(); j++){
				final Membership member=genome.membership(j);
				Counts count=sources.get(member.source);
				if(count==null){
					count=new Counts(member.source);
					sources.put(member.source, count);
					order.add(count);
				}
				count.memberships++;
				memberships++;
				if(count.lastInput!=i){
					count.lastInput=i;
					count.inputs++;
					count.bases=Math.addExact(count.bases, genome.bases);
				}
			}
		}
		if(memberships!=catalog.memberships()){
			throw new IllegalStateException("Per-source membership accounting must conserve the loaded catalog");
		}
		final ByteBuilder b=new ByteBuilder();
		try(Output out=new Output(prefix+SUFFIX)){
			out.writer.println("source\tmembership_rows\tdistinct_input_ids\tnucleotide_bases\tstatus");
			for(final Counts count:order){
				assert(count.inputs>0 && count.inputs<=count.memberships)
					: "Every observed source has members; repeated memberships cannot add input IDs";
				BuildNcrnaReference.field(b.clear(), count.source).tab().append(count.memberships).tab();
				out.writer.print(b.append(count.inputs).tab().append(count.bases)
					.append("\tCATALOG_INPUTS_NOT_INDEPENDENT_ORGANISMS\n"));
			}
		}
	}
	private static final class Counts{
		Counts(final String source_){source=source_;}
		final String source;
		int lastInput=-1, inputs;
		long memberships, bases;
	}
	private NcrnaReferencePanels(){}
	static final String SUFFIX=".panel_population.tsv";
}
