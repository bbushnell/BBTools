package prok;

import java.util.ArrayList;

import map.ObjectIntMap;
import map.ObjectMap;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaReferenceGff.Row;
import prok.SummarizeNcrnaSweep.Table;

/** Optional qualified-support join for search coverage. Startup reads the frozen
 * metadata once and retains only selected inputs and the scored family. After
 * construction all maps/parents are immutable and may be read by file workers.
 * Each parent contributes to every qualified member, never to marginal members.
 * These overlapping counts are not subtype precision or additional union TPs.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaMemberCoverage{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaMemberCoverage(final String membersPath, final String lociPath,
			final String family_, final Binding[] bindings){
		require(bindings!=null && bindings.length>0 && family_!=null, "Member support must bind the actual sweep input population and family");
		family=family_;
		names=names(membersPath, family);
		memberTokens=new byte[names.length][];
		for(int m=0; m<names.length; m++){memberTokens[m]=NcrnaSearchCoverage.token(names[m]);}
		parents=new ArrayList<ObjectMap<String, Parent>>(bindings.length);
		final ObjectIntMap<String> inputs=new ObjectIntMap<String>(String.class);
		for(int i=0; i<bindings.length; i++){
			require(!inputs.contains(bindings[i].id), "Member support cannot duplicate an original input ID");
			inputs.put(bindings[i].id, i);
			parents.add(new ObjectMap<String, Parent>(String.class, Parent.class));
		}
		load(lociPath, inputs);
	}

	/** Member maps use literal UTF-8, exactly as the inventory producer does. */
	static String[] names(final String path, final String family){
		final ArrayList<String> names=new ArrayList<String>();
		final ObjectIntMap<String> seen=new ObjectIntMap<String>(String.class);
		final ObjectIntMap<String> seenModels=new ObjectIntMap<String>(String.class);
		try(Table table=new Table(path)){
			final int f=table.column("family"), m=table.column("member"), model=table.column("reference_model");
			while(table.next()){
				if(!table.literal(f).equals(family)){continue;}
				final String name=table.literal(m), referenceModel=table.literal(model);
				require(!name.isEmpty() && !name.equals(".") && name.indexOf(',')<0 && !seen.contains(name)
					&& !referenceModel.isEmpty() && !referenceModel.equals(".") && !seenModels.contains(referenceModel),
					"Member map requires distinct names and nonmissing, distinct reference models");
				require(names.size()<63, "Coverage uses at most63 positive member-bit positions");
				seen.put(name, names.size());
				seenModels.put(referenceModel, names.size());
				names.add(name);
			}
		}
		require(!names.isEmpty(), "No declared reference members for coverage family: "+family);
		return names.toArray(new String[names.size()]);
	}

	private void load(final String path, final ObjectIntMap<String> inputs){
		final byte[] familyToken=NcrnaSearchCoverage.token(family);
		try(Table table=new Table(path)){
			final int input=table.column("input_id"), id=table.column("parent_id"), f=table.column("family");
			final int seq=table.column("seqid"), start=table.column("start"), stop=table.column("stop"), strand=table.column("strand");
			final int status=table.column("status"), scored=table.column("scored_members"), marginal=table.column("marginal_members");
			byte[] lastInput=null;
			int index=-1;
			while(table.next()){
				if(!table.matches(familyToken, f)){continue;}
				if(lastInput==null || !table.matches(lastInput, input)){
					lastInput=table.parser.parseByteArray(input);
					index=inputs.get(table.text(input));
				}
				if(index<0){continue;}
				final boolean mask=table.parser.termEquals("dont_care", status);
				require(mask || table.parser.termEquals("scored", status), "Coverage support must retain scored/masked inventory status");
				final long qualified=support(table, scored), weak=support(table, marginal);
				require((qualified&weak)==0 && (mask ? qualified==0 : qualified!=0),
					"Masked parents cannot confer qualified coverage, and every scored parent needs qualified support");
				final String parent=table.text(id), contig=table.text(seq), direction=table.text(strand);
				final long from=NcrnaSearchCoverage.count(table, start), to=NcrnaSearchCoverage.count(table, stop);
				require(!parent.isEmpty() && !contig.isEmpty() && from>0 && to>=from && to<=Integer.MAX_VALUE
					&& (direction.equals("+") || direction.equals("-")), "Inventory parent needs exact identity and valid inclusive geometry");
				final ObjectMap<String, Parent> selected=parents.get(index);
				require(!selected.contains(parent), "Repeated inventory parent for one input: "+parent);
				selected.put(parent, new Parent(contig, (int)from, (int)to, direction.charAt(0), mask, qualified));
			}
		}
	}

	/** Compares encoded member tokens in-place; no String/boolean array per support field. */
	private long support(final Table table, final int column){
		if(table.parser.termEquals('.', column)){return 0;}
		table.parser.setBounds(column);
		final byte[] line=table.parser.line();
		final int end=table.parser.b();
		long bits=0;
		for(int begin=table.parser.a(); begin<end;){
			int stop=begin;
			while(stop<end && line[stop]!=','){stop++;}
			int found=-1;
			for(int m=0; m<memberTokens.length; m++){
				final byte[] token=memberTokens[m];
				if(token.length!=stop-begin){continue;}
				boolean same=true;
				for(int i=0; same && i<token.length; i++){same=line[begin+i]==token[i];}
				if(same){
					found=m;
					break;
				}
			}
			require(found>=0 && (bits&(1L<<found))==0, "Inventory support has an unknown, empty or repeated member");
			bits|=1L<<found;
			require(stop==end || stop+1<end, "Trailing empty support member");
			begin=stop+1;
		}
		require(bits!=0, "Missing support must use the explicit dot marker");
		return bits;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Worker Lookup         ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates every reference parent against its inventory twin once per input.
	 * The returned primitive vector is worker-owned, aligned to the live matcher. */
	public long[] bind(final int input, final ArrayList<Row> rows){
		require(input>=0 && input<parents.size() && rows!=null, "Member join needs one declared input and its live reference rows");
		final ObjectMap<String, Parent> selected=parents.get(input);
		require(selected.size()==rows.size(), "Reference and member inventory must contain the same target parents");
		final long[] bits=new long[rows.size()];
		for(int i=0; i<rows.size(); i++){
			final Row row=rows.get(i);
			final Parent parent=selected.get(row.id);
			require(parent!=null && family.equals(row.locus.family) && parent.mask==row.mask
				&& parent.seqid.equals(row.locus.seqid) && parent.start==row.locus.start
				&& parent.stop==row.locus.stop && parent.strand==row.locus.strand,
				"Reference parent status/geometry disagrees with frozen member support: "+row.id);
			bits[i]=parent.bits;
		}
		return bits;
	}

	public int size(){return names.length;}
	public String name(final int member){return names[member];}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}
	private static final class Parent{
		Parent(final String seqid_, final int start_, final int stop_, final char strand_, final boolean mask_, final long bits_){
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
			mask=mask_;
			bits=bits_;
		}
		final String seqid;
		final int start, stop;
		final char strand;
		final boolean mask;
		final long bits;
	}
	private final String family;
	private final String[] names;
	private final byte[][] memberTokens;
	private final ArrayList<ObjectMap<String, Parent>> parents;
}
