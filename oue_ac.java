package th_w;
import org.apache.commons.math3.distribution.BetaDistribution;
//import com.github.davidepattit.XXHash32; // Import XXHash library
import java.util.Random;

import java.util.Random;
import java.util.ArrayList;
public class oue_ac {


    int iter = 100000000;
    int v1 = 1;
    int v2 = 2;
    int binNum = 8;
    double[] epsilon = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 4, 6};


    public oue_ac() throws Exception{
        audit();
    }


    public int getRandomChoice(double[] ue_val, int k) {
        // Sum the elements of ue_val
        int sum = 0;
        for (double val : ue_val) {
            sum += val;
        }

        Random random = new Random();

        // If the sum is 0, return a random integer between 0 and k-1
        if (sum == 0) {
            return random.nextInt(k);
        } else {
            // Create a list to store indices where ue_val[i] == 1
            ArrayList<Integer> indices = new ArrayList<>();
            for (int i = 0; i < ue_val.length; i++) {
                if (ue_val[i] == 1) {
                    indices.add(i);
                }
            }

            // Choose a random index from the collected indices
            return indices.get(random.nextInt(indices.size()));
        }
    }



    public void audit(){
        for(int eps = 0; eps < epsilon.length; eps++) {
            double[] pert1 = new double[binNum];
            double[] pert2 = new double[binNum];
            double[] return1 = new double[binNum];
            double[] return2 = new double[binNum];
            double emp_eps = 0.0;
            int count1 = 0;
            int count2 = 0;
            int test_statistic = 1;
            for (int i = 0; i < iter; i++) {
                return1 = OUE(v1,  epsilon[eps]);
                int predict = getRandomChoice(return1, binNum);
                count1 += (predict == test_statistic) ? 1 : 0;
                return2 = OUE(v2,  epsilon[eps]);
                predict = getRandomChoice(return2, binNum);
                count2 += (predict == test_statistic) ? 1 : 0;
            }

//            for(int i = 0; i < binNum; i++) {
            double[] cp1 = clopperPearson(count1, iter, 0.05);
            double[] cp2 = clopperPearson(count2, iter, 0.05);
            double log_m1 = Math.abs(Math.log(cp1[0] / cp2[1]));
            double log_m2 = Math.abs(Math.log(cp1[1] / cp2[0]));
            emp_eps = Math.abs(Math.max(log_m1, log_m2));
//            }
            System.out.println(count1 + " " + count2 + " " + emp_eps);



        }
    }

    public double[] cal_boundary(int b1, int b2, int binNum){
        double step = (double)(b2 - b1) / binNum;
        double[] boundary = new double[2 * binNum];

        boundary[0] = -1;
        boundary[binNum*2-1] = 1;

        for(int i = 1; i < binNum*2 - 1; i = i + 2){
            boundary[i] = boundary[i - 1] + step;
            boundary[i + 1] = boundary[i];
        }
        return boundary;
    }

    public double laplace(double epsilon) {
        Random random = new Random();
        double sensitivity = 2.0;
        double scale = sensitivity / epsilon;
        // 生成拉普拉斯噪声
        double u = random.nextDouble() - 0.5;
        return -scale * Math.signum(u) * Math.log(1 - 2 * Math.abs(u));
    }

    // 对数据进行扰动
    public double[] perturbation(int value, double epsilon, double[] boundary) {
        double[] encode = new double[binNum];
        // 编码，将特定位置置为 1
        if(value >= boundary[0] && value <= boundary[1])
            encode[0]=1;
        else
            encode[0]=0;
        for(int i = 2; i < binNum * 2; i = i + 2){
            if(value > boundary[i] && value <= boundary[i + 1]){
                encode[i/2] = 1;
            }else{
                encode[i/2] = 0;
            }
        }

        // 添加拉普拉斯噪声
        for (int i = 0; i < binNum; i++) {
            encode[i] = encode[i] + laplace(epsilon);
        }

        return encode;
    }

    // 实现SHE函数
    public double[] SHE(int value, int b1, int b2, double epsilon, double theta) {
        double[] boundary = cal_boundary(b1, b2, binNum);
        double[] perturb_she = perturbation(value, epsilon, boundary);
        return perturb_she;

    }

    // 实现THE函数
    public double[] OUE(int value, double epsilon) {
        double[] ue = new double[binNum];
        for(int i = 0; i < binNum; i++){
            if(value == i)
                ue[i] = 1;
        }
//        double p = Math.exp(epsilon) / (Math.exp(epsilon) + 1);
        double p = 1.0/2;
        double q = 1 / (Math.exp(epsilon) + 1);
        for(int i = 0; i < binNum; i++){
            Random random = new Random();
            double tmp = random.nextDouble();
            if(ue[i] == 1){
                if(tmp > p)
                    ue[i] = 0;
            }else{
                if(tmp < q)
                    ue[i] = 1;
            }
        }
    return ue;
         //support
    }

    public static double[] clopperPearson(int k, int n, double alpha) {
        double lowerBound, upperBound;


        BetaDistribution betaLower = new BetaDistribution(k, n - k + 1);
        lowerBound = betaLower.inverseCumulativeProbability(alpha / 2);


        BetaDistribution betaUpper = new BetaDistribution(k + 1, n - k);
        upperBound = betaUpper.inverseCumulativeProbability(1 - alpha / 2);
        return new double[]{lowerBound, upperBound};
    }
    public static void main(String[] args)  throws Exception {
        new oue_ac();

    }

}

