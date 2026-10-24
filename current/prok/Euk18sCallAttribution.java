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

/** Counts every generic/ribosomal RNA row and binds 18S names to the actual
 * comparison library. This is an attribution audit, not an accuracy grader.
 * @author Brian Bushnell, Raiden */
public final class Euk18sCallAttribution {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, Euk18sCallAttribution.class, false);
		String models=null, calls=null, out=null;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("models")){models=value;}else if(key.equals("calls")){calls=value;}
			else if(key.equals("out")){out=value;}else{throw new IllegalArgumentException("Unknown flag: "+key);}
		}
		require(models!=null && calls!=null && out!=null, "Bind the actual caller library, GFF and audit output");
		Shared.setThreads(1);new Euk18sCallAttribution().run(models, calls, out);Shared.closeStream(pp.outstream);
	}
	public void run(String models, String calls, String output){
		final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
		final ArrayList<String> names=new ArrayList<String>();
		RrnaResourceIO.stream(models, r->{
			require(r.id.startsWith("euk18S_") && r.id.matches("[A-Za-z0-9_]+") && !ids.contains(r.id), "Audit requires unique family-prefixed model names: "+r.id);
			ids.put(r.id, names.size());names.add(r.id);
		});
		require(!names.isEmpty(), "An empty model roster cannot validate emitted 18S calls");
		final long[] modelCounts=new long[names.size()], families=new long[Euk18sM1Grade.FAMILIES.length];
		long rna=0, generic=0, legacy=0;
		final ByteFile in=RrnaResourceIO.open(calls);final LineParser1 p=new LineParser1('\t');
		try{for(byte[] row;(row=in.nextLine())!=null;){
			if(row.length==0 || row[0]=='#'){continue;}
			p.set(row);require(p.terms()==9, "Attribution input must be GFF9");
			final boolean isGeneric=p.termEquals("RNA", 2);
			if(!isGeneric && !p.termEquals("rRNA", 2)){continue;}
			rna++;if(isGeneric){generic++;}
			final String family=RrnaSeedBenchmarkGrade.family(row, p, isGeneric);
			if(family==null){
				final String attrs=p.parseString(8);
				require(!attrs.contains("model:") && (attrs.startsWith("5S,") || attrs.startsWith("16S,") || attrs.startsWith("18S,") || attrs.startsWith("23S,") || attrs.startsWith("28S,")), "Only recognized legacy rRNA may lack a family");
				legacy++;continue;
			}
			int index=-1;for(int i=0; i<families.length; i++){if(Euk18sM1Grade.FAMILIES[i].equals(family)){index=i;break;}}
			require(index>=0, "Unknown emitted RNA family: "+family);families[index]++;
			final String attrs=p.parseString(8);
			final String model=attrs.startsWith("model:") || attrs.contains(",model:") ? RrnaShortlistReplay.model(attrs) : null;
			require(model==null || !model.startsWith("euk18S_") || family.equals("euk18S"), "An euk18S model cannot be attributed to a different family: "+family+" model="+model);
			if(family.equals("euk18S")){
				require(model!=null, "Emitted euk18S call must name its consensus model");
				final int id=ids.get(model);require(id>=0, "Emitted 18S model is absent from its actual library: "+model);modelCounts[id]++;
			}
		}}finally{require(!in.close(), "Read all output rows before asserting attribution conservation");}
		long sum=legacy, euk=0, modelSum=0;
		for(int i=0; i<families.length; i++){sum+=families[i];if(Euk18sM1Grade.FAMILIES[i].equals("euk18S")){euk=families[i];}}
		for(long n:modelCounts){modelSum+=n;}
		require(sum==rna && modelSum==euk, "Family partitions must conserve all RNA rows, and model partitions all euk18S rows");
		RrnaResourceIO.fresh(output, "");final ByteStreamWriter out=RrnaResourceIO.writer(output);final ByteBuilder b=new ByteBuilder();
		try{
			out.println("scope\tname\trows");
			out.print(b.append("TOTAL\tRNA_AND_RRNA\t").append(rna).nl());
			out.print(b.clear().append("TOTAL\tGENERIC_RNA\t").append(generic).nl());
			out.print(b.clear().append("FAMILY\tLEGACY_RRNA\t").append(legacy).nl());
			for(int i=0; i<families.length; i++){out.print(b.clear().append("FAMILY\t").append(Euk18sM1Grade.FAMILIES[i]).tab().append(families[i]).nl());}
			for(int i=0; i<modelCounts.length; i++){out.print(b.clear().append("MODEL\t").append(names.get(i)).tab().append(modelCounts[i]).nl());}
		}finally{RrnaResourceIO.finish(out);}
		System.err.println("EUK18S_ATTRIBUTION_PASS RNA="+rna+" generic="+generic+" euk18S="+euk+" unknown=0");
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}
