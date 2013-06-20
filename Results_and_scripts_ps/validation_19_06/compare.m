
i=8;

err10(1)=abs(kam10(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err10(1)=err10(1)*100;
err10(2)=abs(kam10(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err10(2)=err10(2)*100;

err20(1)=abs(kam20(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err20(1)=err20(1)*100;
err20(2)=abs(kam20(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err20(2)=err20(2)*100;

err30(1)=abs(kam30(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err30(1)=err30(1)*100;
err30(2)=abs(kam30(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err30(2)=err30(2)*100;

err60(1)=abs(kam60(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err60(1)=err60(1)*100;
err60(2)=abs(kam60(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err60(2)=err60(2)*100;

err90(1)=abs(kam90(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err90(1)=err90(1)*100;
err90(2)=abs(kam90(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err90(2)=err90(2)*100;

err120(1)=abs(kam120(i).output_data.Ntot(floor(end/2))-kam200(i).output_data.Ntot(floor(end/2)))/kam200(i).output_data.Ntot(floor(end/2));
err120(1)=err120(1)*100;
err120(2)=abs(kam120(i).output_data.Ntot(floor(end))-kam200(i).output_data.Ntot(floor(end)))/kam200(i).output_data.Ntot(floor(end));
err120(2)=err120(2)*100;

%% Vtot
err10(3)=abs(kam10(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err10(3)=err10(3)*100;
err10(4)=abs(kam10(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err10(4)=err10(4)*100;

err20(3)=abs(kam20(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err20(3)=err20(3)*100;
err20(4)=abs(kam20(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err20(4)=err20(4)*100;

err30(3)=abs(kam30(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err30(3)=err30(3)*100;
err30(4)=abs(kam30(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err30(4)=err30(4)*100;

err60(3)=abs(kam60(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err60(3)=err60(3)*100;
err60(4)=abs(kam60(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err60(4)=err60(4)*100;

err90(3)=abs(kam90(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err90(3)=err90(3)*100;
err90(4)=abs(kam90(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err90(4)=err90(4)*100;

err120(3)=abs(kam120(i).output_data.Vtot(floor(end/2))-kam200(i).output_data.Vtot(floor(end/2)))/kam200(i).output_data.Vtot(floor(end/2));
err120(3)=err120(3)*100;
err120(4)=abs(kam120(i).output_data.Vtot(floor(end))-kam200(i).output_data.Vtot(floor(end)))/kam200(i).output_data.Vtot(floor(end));
err120(4)=err120(4)*100;
