package test.prok;

import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;

import java.util.ArrayList;
import java.util.Locale;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import shared.Timer;

/** Experimental full CallGenes driver; only explicit mitochondrial resources
 * enable the new refinement. Every ordinary family flag is forwarded unchanged.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointCallGenes {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, MitoEndpointCallGenes.class, false);
		try{run(pp.args);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Parses only experiment options; the ordinary caller owns all other flags. */
	private static void run(final String[] args){
		final ArrayList<String> forwarded=new ArrayList<String>();
		String net=null, tables=null, pin=null;
		int radius=192;
		boolean radiusSet=false;
		for(String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("mitonet")){net=b;}
			else if(a.equals("mitotables")){tables=b;}
			else if(a.equals("mitolibrarysha80")){pin=b;}
			else if(a.equals("mitoradius")){
				radius=Integer.parseInt(b);
				radiusSet=true;
			}else{forwarded.add(arg);}
		}
		final boolean enabled=net!=null;
		require(enabled==(tables!=null) && enabled==(pin!=null) && (enabled || !radiusSet),"Mito refinement needs explicit net/tables/library pin together; omitted resources mean exact NN-off control");
		final Timer timer=new Timer();
		final CallGenes caller=new CallGenes(forwarded.toArray(new String[0]));
		final String netPath=net, tablePath=tables, libraryPin=pin;
		final int searchRadius=radius;
		caller.run(timer, enabled ? () -> configure(netPath, tablePath, libraryPin, searchRadius) : null);
		System.err.println("MITO_ENDPOINT_CALLGENES_PASS enabled="+enabled+" radius="+(enabled ? radius : 0));
	}

	/** Called inside CallGenes' owned scope before any calling worker starts. */
	private static void configure(final String net, final String tables, final String pin, final int radius){
		NcrnaFamily selected=null;
		for(NcrnaFamily family : GeneCaller.ncrnaFamilies){
			if("16S".equals(family.name)){
				require(selected==null, "Only one 16S family may own the mitochondrial endpoint resources");
				selected=family;
			}
		}
		require(selected!=null, "Enable 16S before supplying mitochondrial endpoint resources");
		selected.boundaryRefinerFactory=MitoEndpointCaller.Resources.load(selected, net, tables, pin, radius);
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
