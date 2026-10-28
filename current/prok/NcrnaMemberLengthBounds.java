package prok;

import java.util.ArrayList;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Exports the existing effective length_bounds contract from native admitted
 * construction-member extrema: floor(0.8*minimum), ceil(1.2*maximum).
 * The immutable input report retains raw extrema and membership counts.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaMemberLengthBounds {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaMemberLengthBounds.class, false);
		String input=null, refs=null, output=null, family=null;int expected=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){input=v;}else if(k.equals("ref")){refs=v;}else if(k.equals("out")){output=v;}else if(k.equals("family")){family=v;}else if(k.equals("expected")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}}
		require(input!=null && refs!=null && output!=null && family!=null && family.matches("RF[0-9]{5}") && expected>0, "Bind native membership report, exact consensus library, family and model count");Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{run(input, refs, output, family, expected);}finally{Shared.closeStream(pp.outstream);}
	}
	private static void run(String input, String refs, String output, String family, int expected){
		final ArrayList<Model> models=new ArrayList<Model>();final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
		RrnaResourceIO.stream(refs, r->{final String id=RrnaResourceIO.first(r.id);require(id.startsWith("cobalamin_"+family+'_') && !ids.contains(id) && r.bases!=null && r.length()>0, "Consensus names must uniquely belong to the declared family");ids.put(id, models.size());models.add(new Model(r.id, r.length()));});
		require(models.size()==expected, "Every final consensus needs a bound row");
		final ByteFile in=RrnaResourceIO.open(input);final LineParser1 p=new LineParser1('\t');int rows=0;long members=0;
		try{RrnaResourceIO.header(in, HEADER, input);for(byte[] row;(row=in.nextLine())!=null;){p.set(row);require(p.terms()==19 && ids.contains(p.parseString(0)) && p.termEquals(family, 1), "Membership report must match the exact final model set and emitted family");
			final Model m=models.get(ids.get(p.parseString(0)));final int n=p.parseInt(3), lo=p.parseInt(4), hi=p.parseInt(8);
			require(!m.seen && p.parseInt(2)==m.length && n>0 && lo>0 && hi>=lo && (p.termEquals("PURE", 10) || p.termEquals("MIXED", 10)), "Empty, duplicate or changed models cannot get invented length limits");
			int previous=lo;for(int col=5; col<=8; col++){final int value=p.parseInt(col);require(value>=previous && value<=hi, "Member length order statistics must be monotonic");previous=value;}
			long labels=0;for(int col=11; col<=17; col+=2){final int count=p.parseInt(col);require(count>=0, "Biological membership counts cannot be negative");labels+=count;}
			require(labels==n, "Biological labels must conserve construction memberships");m.seen=true;
			m.min=NcrnaMitoLengthProfile.lowerBound(lo);m.max=NcrnaMitoLengthProfile.upperBound(hi);rows++;members+=n;
		}}finally{require(!in.close(), "Membership report read failed");}
		require(rows==expected, "Missing model bounds cannot silently become defaults");for(Model m:models){require(m.seen, "Missing model: "+m.id);}
		RrnaResourceIO.fresh(output, "");final ByteStreamWriter out=RrnaResourceIO.writer(output);final ByteBuilder b=new ByteBuilder();
		try{out.println("model\tmodel_length\tmin_length\tmax_length");for(Model m:models){out.print(b.clear().append(m.id).tab().append(m.length).tab().append(m.min).tab().append(m.max).nl());}}finally{RrnaResourceIO.finish(out);}
		System.err.println("NCRNA_MEMBER_LENGTH_BOUNDS_PASS models="+models.size()+" memberships="+members+"; existing outward floor(0.8*rawMin)/ceil(1.2*rawMax), raw report retained; never padded-window length");
	}
	private static final class Model {Model(String id_, int length_){id=id_;length=length_;}final String id;final int length;int min,max;boolean seen;}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String HEADER="model_id\temitted_family\tconsensus_length\tconstruction_members\tmember_min\tmember_p10\tmember_median\tmember_p90\tmember_max\temitted_family_fraction\tpurity\tRF00174_members\tRF00174_fraction\tRF01482_members\tRF01482_fraction\tRF01689_members\tRF01689_fraction\tB12clade_members\tB12clade_fraction";
}
