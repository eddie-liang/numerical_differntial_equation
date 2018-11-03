function [ ] = output(tm, er, dt, t1, u1, U1)

  disp([' ']);
  disp(['Computed solution in ', num2str(tm), ' secs'])
  disp(['With dt = ', num2str(dt), ' max error is ', num2str(er)])
  disp(' ');

  % Plot the solution
  figure(1); clf;
  plot(t1,u1,'-r','LineWidth',4);
  hold on;
  plot(t1,U1,'-o','MarkerSize',10);
  hold off;
  legend('Exact','Euler');
  xlabel('time');
  title('Exp Solution');
  drawnow;

  figure(2); clf;
  plot(t1, abs(u1 - U1),'-ob','LineWidth',4,'MarkerSize',10)
  xlabel('time');
  title('Error in Exp Solution');
  drawnow;

end
