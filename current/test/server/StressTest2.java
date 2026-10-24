package test.server;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.time.Duration;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import parse.Parse;
import shared.Timer;

/**
 * Bounded fixed-pool Java 11 HttpClient load helper for finite local endpoints.
 * Uses 10 second connect and 30 second request timeouts so a silent fixture
 * endpoint cannot hang the test indefinitely.
 */
public class StressTest2 {

	private static final Duration CONNECT_TIMEOUT=Duration.ofSeconds(10);
	private static final Duration REQUEST_TIMEOUT=Duration.ofSeconds(30);

	/**
	 * Sends blocking GET requests to an explicit endpoint from a bounded worker
	 * pool.
	 * @param args endpoint, optional iteration count, optional max concurrency
	 */
	public static void main(String[] args){
		Config config=new Config(args);
		ExecutorService executor=Executors.newFixedThreadPool(config.maxConcurrent);
		ExecutorCompletionService<Integer> service=new ExecutorCompletionService<Integer>(executor);
		HttpClient client=HttpClient.newBuilder()
				.version(HttpClient.Version.HTTP_2)
				.connectTimeout(CONNECT_TIMEOUT)
				.build();
		HttpRequest request=HttpRequest.newBuilder()
				.uri(config.endpoint)
				.timeout(REQUEST_TIMEOUT)
				.build();
		Stats stats=new Stats();
		Timer t=new Timer();

		try{
			long submitted=0;
			for(int i=0; i<config.maxConcurrent && submitted<config.iterations; i++, submitted++){
				submit(service, client, request);
			}
			for(long received=0; received<config.iterations; received++){
				try{
					stats.record(service.take().get().intValue());
				}catch(ExecutionException e){
					stats.exceptions++;
				}
				if(submitted<config.iterations){
					submit(service, client, request);
					submitted++;
				}
			}
		}catch(InterruptedException e){
			Thread.currentThread().interrupt();
			stats.exceptions++;
		}finally{
			executor.shutdownNow();
		}

		t.stop("Completed "+stats.completed+" of "+config.iterations+" pooled queries in ");
		stats.finish(config.iterations);
	}

	/** Submits one blocking request. */
	private static void submit(ExecutorCompletionService<Integer> service, final HttpClient client,
			final HttpRequest request){
		service.submit(new Callable<Integer>() {
			@Override
			public Integer call() throws Exception {
				return Integer.valueOf(client.send(request, HttpResponse.BodyHandlers.discarding()).statusCode());
			}
		});
	}

	/** Parsed command-line options. */
	private static class Config {
		Config(String[] args){
			if(args.length<1 || args.length>3){throw new IllegalArgumentException(
					"Usage: java test.server.StressTest2 <endpoint> [iterations=10000] [maxConcurrent=availableProcessors]");}
			endpoint=URI.create(args[0]);
			iterations=(args.length>1 ? Parse.parseKMG(args[1]) : 10000);
			maxConcurrent=(args.length>2 ? parsePositiveInt(args[2], "maxConcurrent") :
				Math.max(1, Runtime.getRuntime().availableProcessors()));
			if(iterations<1){throw new IllegalArgumentException("Iterations must be positive: "+iterations);}
			if(maxConcurrent<1){throw new IllegalArgumentException("maxConcurrent must be positive: "+maxConcurrent);}
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
		final int maxConcurrent;
	}

	/** Response counters and failure policy independent of Java assertions. */
	private static class Stats {
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
