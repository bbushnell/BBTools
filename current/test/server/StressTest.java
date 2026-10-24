package test.server;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.time.Duration;
import java.util.ArrayDeque;
import java.util.concurrent.CancellationException;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;

import parse.Parse;
import shared.Timer;

/**
 * Bounded async Java 11 HttpClient load helper for finite local endpoints.
 * Uses 10 second connect and 30 second request timeouts so a silent fixture
 * endpoint cannot hang the test indefinitely.
 */
public class StressTest {

	private static final Duration CONNECT_TIMEOUT=Duration.ofSeconds(10);
	private static final Duration REQUEST_TIMEOUT=Duration.ofSeconds(30);

	/**
	 * Sends asynchronous GET requests to an explicit endpoint while limiting the
	 * number of outstanding requests.
	 * @param args endpoint, optional iteration count, optional max outstanding
	 */
	public static void main(String[] args){
		Config config=new Config(args, 64);
		HttpClient client=HttpClient.newBuilder()
				.version(HttpClient.Version.HTTP_2)
				.connectTimeout(CONNECT_TIMEOUT)
				.build();
		HttpRequest request=HttpRequest.newBuilder()
				.uri(config.endpoint)
				.timeout(REQUEST_TIMEOUT)
				.build();
		ArrayDeque<CompletableFuture<Integer>> pending=new ArrayDeque<CompletableFuture<Integer>>();
		Stats stats=new Stats();
		Timer t=new Timer();

		for(long i=0; i<config.iterations; i++){
			pending.addLast(client.sendAsync(request, HttpResponse.BodyHandlers.discarding())
					.thenApply(response -> Integer.valueOf(response.statusCode())));
			if(pending.size()>=config.maxOutstanding){stats.record(pending.removeFirst());}
		}
		while(!pending.isEmpty()){stats.record(pending.removeFirst());}

		t.stop("Completed "+stats.completed+" of "+config.iterations+" async queries in ");
		stats.finish(config.iterations);
	}

	/** Parsed command-line options. */
	private static class Config {
		Config(String[] args, int defaultOutstanding){
			if(args.length<1 || args.length>3){throw new IllegalArgumentException(
					"Usage: java test.server.StressTest <endpoint> [iterations=10000] [maxOutstanding=64]");}
			endpoint=URI.create(args[0]);
			iterations=(args.length>1 ? Parse.parseKMG(args[1]) : 10000);
			maxOutstanding=(args.length>2 ? parsePositiveInt(args[2], "maxOutstanding") : defaultOutstanding);
			if(iterations<1){throw new IllegalArgumentException("Iterations must be positive: "+iterations);}
			if(maxOutstanding<1){throw new IllegalArgumentException("maxOutstanding must be positive: "+maxOutstanding);}
		}
		private static int parsePositiveInt(String value, String name){
			long x=Parse.parseKMG(value);
			if(x<1 || x>Integer.MAX_VALUE){
				throw new IllegalArgumentException(name+" must be in [1, "+Integer.MAX_VALUE+"]: "+x);
			}
			return (int)x;
		}
		final URI endpoint;
		final long iterations;
		final int maxOutstanding;
	}

	/** Response counters and failure policy independent of Java assertions. */
	private static class Stats {
		void record(CompletableFuture<Integer> future){
			try{
				record(future.join().intValue());
			}catch(CompletionException e){
				exceptions++;
			}catch(CancellationException e){
				exceptions++;
			}
		}
		void record(int status){
			completed++;
			if(status<200 || status>299){non2xx++;}
		}
		void finish(long expected){
			System.out.println("completed\t"+completed);
			System.out.println("non2xx\t"+non2xx);
			System.out.println("exceptions\t"+exceptions);
			if(completed!=expected || non2xx>0 || exceptions>0){
				throw new IllegalStateException("completed="+completed+", expected="+expected+
						", non2xx="+non2xx+", exceptions="+exceptions);
			}
		}
		long completed;
		long non2xx;
		long exceptions;
	}
}
